#include <string.h> 		/* memset */
#include <stdlib.h>
#include <argp.h> /*parsing arguments*/
#include <TTree.h>
#include <TFile.h>

const int C_ROC_READ_LIMIT = 5000000;/*for debugging - reading only few RO cycles*/
const int C_MAX_BIF_EVENTS = 10*1000*1000;

static char doc[] = "Statistics printout tool for RAW data dumped from EUDAQ BIF producer.  For more info try --help";
static char args_doc[] = "";

static struct argp_option options[] =
{
    { "spiroc_raw_file", 'w', "SPIROC_RAW_FILE", 0, "raw file for the statistics data from AHCAL RAW file saved by EUDAQ" },
    { "bif_raw_file", 'b', "BIF_FILE", 0, "filename/path of the BIF raw data file" },
    { "dwc_root_file", 'd', "DWC_ROOT_FILE", 0, "filename/path of the DWC root data file" },
    
    { 0 } };

/* Used by main to communicate with parse_opt. */
struct arguments_t {
    char *bif_filename;
    char *spiroc_raw_filename;
    char *dwc_root_filename;
};
struct arguments_t arguments;

void arguments_init(struct arguments_t* arguments) {
    /* Default values. */
    arguments->spiroc_raw_filename = NULL;
    arguments->bif_filename = NULL;
    arguments->dwc_root_filename = NULL;
    //   arguments->print_bif_start_phases = 0;
}

void arguments_print(struct arguments_t* arguments) {
    printf("#BIF_data_file=\"%s\"\n", arguments->bif_filename);
    printf("#SPIROC_RAW_data_file=\"%s\"\n", arguments->spiroc_raw_filename);
    printf("#DWC_ROOT_data_file=\"%s\"\n", arguments->dwc_root_filename);
}

/* Parse a single option. */
static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
    struct arguments_t *arguments = (struct arguments_t *)state->input;
    
    switch (key) {
        case 'b':
            arguments->bif_filename = arg;
            break;
        case 'w':
            arguments->spiroc_raw_filename = arg;
            break;
        case 'd':
            arguments->dwc_root_filename = arg;
            break;
        case ARGP_KEY_END:
            if (arguments->bif_filename == NULL) {
                argp_error(state, "missing BIF data filename!\n");
            }
            if (arguments->spiroc_raw_filename == NULL) {
                argp_error(state, "missing AHCAL RAW filename!\n");
            }
            if (arguments->dwc_root_filename == NULL) {
                argp_error(state, "missing DWC ROOT filename!\n");
            }
            //			argp_usage(state);
            //			argp_state_help(state, state->err_stream, ARGP_HELP_USAGE);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

typedef struct {
    u_int64_t tdc; //48-bit+5bit finestamp
    u_int32_t trig_count; //bif trigger counter
    u_int32_t ro_cycle; //derived from the Readout cycle corresponding to the start of acquisition
} BIF_record_t;

int update_counter_modulo(unsigned int oldvalue, unsigned int newvalue_modulo, unsigned int modulo, unsigned int max_backwards){
    unsigned int newvalue = oldvalue - max_backwards;
    unsigned int mask = modulo - 1;
    if ((newvalue & mask) > (newvalue_modulo & mask)) {
        newvalue += modulo;
    }
    newvalue = (newvalue & (~mask)) | (newvalue_modulo & mask);
    return newvalue;
}


int load_timestamps_from_ahcal_raw(FILE * fp, BIF_record_t * ahcal_data) {

    u_int32_t ROC = ahcal_data->ro_cycle;
    u_int64_t TS = 0;
    
    unsigned char minibuf[8];
    while (1) {
        if (fread(minibuf, 1, 1, fp) <= 0) return -1;
        if (minibuf[0] != 0xCD) continue;
        if (fread(minibuf, 1, 1, fp) <= 0) return -1;
        if (minibuf[0] != 0xCD) continue;
        if (fread(minibuf, sizeof(minibuf), 1, fp) <= 0) return -1;
        int length = (int) minibuf[0] + ((int) minibuf[1] << 8);
        if ((minibuf[0] != 0x10) || (minibuf[1] != 0x0) || (minibuf[7] != 0x08)) {
            //not packets we want
            fseek(fp, length, SEEK_CUR);
            continue;
        }
        int newROC = (int) minibuf[2];
        if (fread(minibuf, sizeof(minibuf), 1, fp) <= 0) return -1;
        if ((minibuf[0] != 0x45) || (minibuf[1] != 0x4D) || (minibuf[2] != 0x49) || (minibuf[3] != 0x54)) {
            fseek(fp, 8, SEEK_CUR);
            continue;
        }
        int type = minibuf[4];
        
        int trigid = ((int) minibuf[6]) + (((int) minibuf[7]) << 8); //possible trigid

        if (fread(minibuf, sizeof(minibuf), 1, fp) <= 0) return -1;
        if ((minibuf[6] != 0xAB) || (minibuf[7] != 0xAB)) continue;
        TS = (u_int64_t) minibuf[0] +
        ((u_int64_t) minibuf[1] << 8) +
        ((u_int64_t) minibuf[2] << 16) +
        ((u_int64_t) minibuf[3] << 24) +
        ((u_int64_t) minibuf[4] << 32) +
        ((u_int64_t) minibuf[5] << 40);
        
        if (type != 0x10) {
            if (type == 0x01) { //start acq
                ROC = update_counter_modulo(ROC, newROC, 256, 10);
            } else if (type == 0x02) { //stop acq
            }
            continue;
        }
        int increment = (newROC - ROC) & 0xFF;
        if (increment > 50) {
            fprintf(stdout, "#ERROR wrong increment of ROC: %d\n", increment);
            //continue;
        }
        ROC = ROC + increment;
        
        ahcal_data->ro_cycle = ROC;
        ahcal_data->tdc = TS;
        ahcal_data->trig_count = trigid;
        return 1;
    }

}

int load_bif_data(FILE * fp, BIF_record_t * bif_data, BIF_record_t * bif_first_data) {
    u_int32_t time_h = 0, time_l = 0;
    u_int64_t time = 0;
    u_int64_t finetime_trig = 0;
    u_int32_t trig_counter = 0;
    u_int8_t trig_details[4];
    u_int64_t shutter_cnt = bif_data->ro_cycle;
    u_int16_t details = shutter_cnt & 0x0FFF;
    u_int16_t details_last = 0;
    unsigned char minibuf[8];
    
    while (-1) {
        if (fread(minibuf, sizeof(minibuf), 1, fp) <= 0)
            return -1;
        time_h = *((u_int32_t *) minibuf + 1) & 0x0000FFFF;
        time_l = *(u_int32_t *) minibuf;
        time = ((u_int64_t) time_h << 32) | time_l;
        details_last = details;
        details = (*((u_int16_t *) minibuf + 3)) & 0x0FFF;
        
        switch (minibuf[7] >> 4) {
            case 0:
            case 1:
                if (fread(minibuf, sizeof(minibuf), 1, fp) <= 0)                    return -1;
                trig_counter = *((u_int32_t *) minibuf + 0);
                if (bif_first_data->trig_count == 0) {
                    bif_first_data->trig_count = trig_counter;
                }
                trig_details[0] = minibuf[4];
                finetime_trig = (time << 5) | ((trig_details[0] + 0x18) & 0x1F);
                
                bif_data->ro_cycle = shutter_cnt;
                bif_data->tdc = finetime_trig;
                bif_data->trig_count = trig_counter;
                return 1;
            case 2:
                break;
            case 3:
                if ((details_last == 4095) && (details == 0)) {
                    /*when the shutter counter overflows, we have to increment the shutter counter properlyew*/
                    shutter_cnt += 4096;
                }
                shutter_cnt &= 0xFFFFFFFFFFFFF000;
                shutter_cnt |= details;
                if (bif_first_data->ro_cycle == 0) {
                    bif_first_data->ro_cycle = shutter_cnt;
                }
                break;
            
            default:
                fprintf(stderr, "unknown BIF packet\n");
                return 0;
                break;
        }
    }
}


int main(int argc, char **argv) {
    /* Default values. */
    arguments_init(&arguments);
    BIF_record_t bif_data = (BIF_record_t) {0, 0, 0};
    BIF_record_t bif_first_data = (BIF_record_t) {0, 0, 0};
    BIF_record_t ahcal_data = (BIF_record_t) {0, 0, 0};
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    arguments_print(&arguments);

    int i;
    u_int64_t TS_diff = 0;
    
    FILE *bif_file;
    if (!(bif_file = fopen(arguments.bif_filename, "r"))) {
        perror("Unable to open the bif file\n");
        return -1;
    }
    if (load_bif_data(bif_file, &bif_data, &bif_first_data) != 1) {
        perror("Unable to read the bif file\n");
        fclose(bif_file);
        return -1;
    }

    FILE *ahcal_file;
    if (!(ahcal_file = fopen(arguments.spiroc_raw_filename, "r"))) {
        perror("Unable to open the ahcal file\n");
        fclose(bif_file);
        return -1;
    }
    if (load_timestamps_from_ahcal_raw(ahcal_file, &ahcal_data) != 1) {
        perror("Unable to read the ahcal file\n");
        fclose(bif_file);
        fclose(ahcal_file);
        return -1;
    }
    
    TFile *dwc_file = new TFile(arguments.dwc_root_filename,"read");
    TTree *dwc_tree = (TTree*) dwc_file->Get("DelayWireChambers");
    
    u_int32_t dwc_event;
    Long64_t dwc_time;
    dwc_tree->SetBranchAddress("event",&dwc_event);
    dwc_tree->SetBranchAddress("timeSinceStart",&dwc_time);
    
    int nEvents = dwc_tree->GetEntries();
    int dwc_index=0;
    
    int event_num = 0, matched = 0;
    int match_index = 1;
    
//    FILE *output = fopen("combined.txt","w");
    TFile *output = new TFile("combined.root","recreate");
    TTree *outtree = new TTree("combined","combined events");
    int ROC;
    u_int32_t bif_trig, ahc_trig;
    u_int64_t biftime;

//    fprintf(output,"ROC\tBifTrg#\tAHCTrg#\ttdc (in BIF)\n");
    outtree->Branch("ROC",&ROC,"ROC/I");
    outtree->Branch("bif_Trig",&bif_trig,"bif_Trig/i");
    outtree->Branch("ahc_Trig",&ahc_trig,"ahc_Trig/i");
    outtree->Branch("bif_Time",&biftime,"bif_Time/l");
    outtree->Branch("dwc_Trig",&dwc_event,"dwc_Trig/i");
    outtree->Branch("dwc_Time",&dwc_time,"dwc_Time/L");

    while (-1) {
        biftime = bif_data.tdc>>5;
        if (TS_diff == 0) TS_diff = biftime - ahcal_data.tdc;
        dwc_tree->GetEntry(dwc_index++);

        if (biftime - ahcal_data.tdc < TS_diff-1) {
            //bif data only
            //fprintf(output,"%d\t%d\t-------\t%ld\n",bif_data.ro_cycle-bif_first_data.ro_cycle+1, bif_data.trig_count-bif_first_data.trig_count, biftime);
            ROC = bif_data.ro_cycle-bif_first_data.ro_cycle+1;
            bif_trig = bif_data.trig_count-bif_first_data.trig_count;
            ahc_trig = -1;
            outtree->Fill();
            event_num++;
            if (load_bif_data(bif_file, &bif_data, &bif_first_data) != 1) {
                break;
            }
        } else if (biftime - ahcal_data.tdc > TS_diff+1) {
            //ahcal data only
            //fprintf(output,"%d\t-------\t%d\t%ld\n",ahcal_data.ro_cycle, ahcal_data.trig_count, ahcal_data.tdc+TS_diff);
            ROC = ahcal_data.ro_cycle;
            bif_trig = -1;
            ahc_trig = ahcal_data.trig_count;
            biftime = ahcal_data.tdc+TS_diff;
            outtree->Fill();
            event_num++;
            if (load_timestamps_from_ahcal_raw(ahcal_file, &ahcal_data) != 1) {
                break;
            }
        } else {
            //matched data
            //fprintf(output,"%d\t%d\t%d\t%ld\n",bif_data.ro_cycle-bif_first_data.ro_cycle+1, bif_data.trig_count-bif_first_data.trig_count, ahcal_data.trig_count, biftime);
            ROC = bif_data.ro_cycle-bif_first_data.ro_cycle+1;
            bif_trig = bif_data.trig_count-bif_first_data.trig_count;
            ahc_trig = ahcal_data.trig_count;
            outtree->Fill();
            event_num++;
            matched++;
            if (load_bif_data(bif_file, &bif_data, &bif_first_data) != 1) {
                break;
            }
            if (load_timestamps_from_ahcal_raw(ahcal_file, &ahcal_data) != 1) {
                break;
            }
        }
        if (event_num == 20 && matched <= 10) {
            match_index++;
            perror("event matching failed\n");
            printf("new matching index: %d\n",match_index);
            //fprintf(output,"\nmatching failed\nmatching index set to %d\n",++match_index);
            fseek(bif_file, 0, SEEK_SET);
            fseek(ahcal_file, 0, SEEK_SET);
            load_bif_data(bif_file, &bif_data, &bif_first_data);
            load_timestamps_from_ahcal_raw(ahcal_file, &ahcal_data);
            for (i=0; i<match_index/2; i++) {
                if (match_index%2) load_bif_data(bif_file, &bif_data, &bif_first_data);
                else load_timestamps_from_ahcal_raw(ahcal_file, &ahcal_data);
            }
            outtree->Reset();
            event_num=0;
            matched=0;
            TS_diff=0;
            dwc_index=0;
        }
        if (dwc_index==nEvents) break;
    }
    
    outtree->Write();
    output->Close();

    dwc_file->Close();
    
    fclose(bif_file);
    fclose(ahcal_file);
    //fclose(output);
				
    exit(0);
}
