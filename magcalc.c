/* Program: galamags
   Usage: Calculate observed AB magnitudes.
   Author: Chuanwu Liu.
   Instit: The University of Melbourne.
   Date: 08/08/2014
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define APPARENT 1 // Calculate apparent magnitudes? =1: Yes. =0: No.
#define ABSOLUTE 1 // Calculate absolute magnitudes? =1: Yes. =0: No.
#define LYALPHA 1 // Employ Lyman-alpha absorption for apparent magnitudes? =1: Yes. =0: No. 
#define DUSTEXTINC 1 // Dust extinction;

#define SNAP_MAX 97
#define NSMOOTH 24
#define NWAVES 1221
#define FILTERLOW 1700 // Include FILTERLOW
#define FILTERUP 17000 // Do not include FILTERUP

const float Omega_M = 0.308;
const float Omega_k = 0.0;
const float Omega_L = 0.692;
const float h0 = 0.678;
const float Tdelay = 0E6;

char *fltlist[9][2] = {{"V606", "/home/cliu/Works/filters/f606w.UVIS1.tab"}, {"i775", "/home/cliu/Works/filters/f775w.UVIS1.tab"},
                       {"I814", "/home/cliu/Works/filters/f814w.UVIS1.tab"}, {"z850", "/home/cliu/Works/filters/f850lp.UVIS1.tab"}, 
                       {"Y098", "/home/cliu/Works/filters/f098m.IR.tab"}, {"Y105", "/home/cliu/Works/filters/f105w.IR.tab"},
                       {"J125", "/home/cliu/Works/filters/f125w.IR.tab"}, {"H160", "/home/cliu/Works/filters/f160w.IR.tab"},
                       {"[3.6]", "/home/cliu/Works/filters/filter-0198.asc"}};
int bands[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

char *rest_title[2] = {"M1500", "M1600"};
float rest_bands[2][2] = {{1500.0, 100.0}, {1600.0, 100.0}}; // rest band center and width.

char *model_seds[5]= {"/home/cliu/Works/spectra/0.05solarmetal/0.05solarmetal.spectrum1",
                    "/home/cliu/Works/spectra/Geneva_std_0.004/Geneva_std_0.004.spectrum1",
                    "/home/cliu/Works/spectra/Geneva_std_0.008/Geneva_std_0.008.spectrum1",
                    "/home/cliu/Works/spectra/Geneva_std_0.020/Geneva_std_0.020.spectrum1",
                    "/home/cliu/Works/spectra/Geneva_std_0.040/Geneva_std_0.040.spectrum1"};

float metal_list[5] = {0.001, 0.004, 0.008, 0.020, 0.040};
float metal_rank[5] = {0, 3, 7, 19, 39};

double Energy[NWAVES];
double Flux[NWAVES];
double spectra[40][NSMOOTH*(SNAP_MAX - 1) + 1][NWAVES];
double (*spectrum)[NWAVES];
float Wave[NWAVES];
float Waveobs[NWAVES];
float Dust[NWAVES], Absorb[NWAVES];
double Flxobs[NWAVES];
double Thrutot[9];
double Thrutot_rest[2];
float Thru[9][FILTERUP - FILTERLOW];
float Thru_rest[2][FILTERUP - FILTERLOW];
float z[SNAP_MAX], lbtime[SNAP_MAX];
float z0, metal;
int Nsnaps, isnap, Nprog, Zrank, iZ, Snap_stop, Snap_start;
int kmetal[SNAP_MAX], imetal;
int main(void)
{
    FILE *openfile(char *fname, char *mode);
    double read_thru(char *fltname, float *thru);
    double read_thru36(char *fltname, float *thru);
    float read_data(FILE *fp, float *data);
    double creat_thru(float *thru, float bcenter, float bwidth);
    double d_lumin(float z);
  //void dust_curve(float *wave);
    void absorb_table(float z, float *waveobs);
    void thru_filter(float *thru, float *wave, double *flx, double *flxobs);
    void thru_filter02(float *thru, float *wave, double *flx, double *flxobs);
    void thru_filter36(float *thru, float *wave, double *flx, double *flxobs);
    float abmag(float *wave, double *flux, double thoutot);

    FILE *fp, *fp0, *fp1;
    time_t tmstart, tmend;
  //double spectrum[NSMOOTH*(SNAP_MAX - 1) + 1][NWAVES];
    float spectexp[NWAVES];
    float newmass[SNAP_MAX], metallicity[SNAP_MAX], newstars;
    double dl;
    float tb, zi, zplus, dz, tmp, starage, sfr, colorexcess;
    register int i, j, k, l, n;


    fputs("================ START =================", stdout);
    tmstart = time(NULL);
    printf("Run Start: %s\n", asctime(localtime(&tmstart)));

    fputs("====== Read the data of galaxies ======", stdout);

    fp0 = fopen("sfrhist.bin", "rb");
    // Read the number of snapshots, z list and look back time list.
    fread(&Snap_start, 4, 1, fp0);
    fread(&Snap_stop, 4, 1, fp0);
    Nsnaps = Snap_stop - Snap_start;
    fread(z, 4, Snap_stop, fp0);
    fread(lbtime, 4, Snap_stop, fp0);
    z0 = z[Snap_stop - 1];
        printf("Snap_start = %d, Snap_stop = %d\n", Snap_start, Snap_stop);
    for (i = 0; i < Snap_stop; i++) {
        printf("z = %.3f  lbt = %.3f\n", z[i], lbtime[i]);
    }


    fputs("=========== Read the wavelength ==========\n", stdout);
    fp = openfile(model_seds[0], "r");
    for (i = 0; i < 6; ++i)
        while(fgetc(fp) != '\n')
            ;

    for (k = 0; k < NWAVES; ++k) {
        fseek(fp, 16L, 1);
        fscanf(fp, "%f%*[^\n]", Wave + k);
    }
    fclose(fp);

    fputs("======= Read the spectra at different look back time =======\n", stdout);

    for (iZ = 0; iZ < 5; iZ ++) {

        Zrank = metal_rank[iZ];
     
        fp = openfile(model_seds[iZ], "r");

        fseek(fp, 0L, 0);
        for (i = 0; i < 6; ++i)
            while(fgetc(fp) != '\n')
                ;

        tb = read_data(fp, spectexp);
        for (i = Nsnaps - 1; i > 0; --i) {
            printf("Info: Processing Z %.3f Snap %d ...\n", metal_list[iZ], i);
            for (j = 0; j < NSMOOTH; ++j) {
                starage = 1e6*(lbtime[i + Snap_start] + (lbtime[i + Snap_start -1] - lbtime[i + Snap_start])*j/NSMOOTH -
                                                                                        lbtime[Snap_stop - 1]) + Tdelay;
                while (1) {
                    if (tb >= starage && tb < starage + 1.4E6) {
                        for (k = 0; k < NWAVES; ++k)
                            spectra[Zrank][i*NSMOOTH - j][k] = pow(10.0, spectexp[k]);
                        break;
                    }
                    else
                        tb = read_data(fp, spectexp);

                }
            }
        }

        printf("Info: Processing Snap %d ...\n", 0);
        starage = 1E6*(lbtime[Snap_start] - lbtime[Snap_stop - 1]) + Tdelay;
        while (1) {
            if (tb >= starage && tb < starage + 1.4E6) {
                for (k = 0; k < NWAVES; ++k)
                    spectra[Zrank][0][k] = pow(10.0, spectexp[k]);
                break;
            }
            else
                tb = read_data(fp, spectexp);
        }
        fclose(fp);
    }

    // evaluate the spectrum using linear interpolation.
    for (iZ = 0; iZ < 40; iZ++) {
        if (iZ > 0 && iZ < 3)
             for (j = 0; j < NSMOOTH*(SNAP_MAX - 1) + 1; ++j)
                  for (k = 0; k < NWAVES; ++k)
                       spectra[iZ][j][k] = spectra[0][j][k] + (spectra[3][j][k] - spectra[0][j][k])/3*iZ;
        else if (iZ > 3 && iZ < 7)
             for (j = 0; j < NSMOOTH*(SNAP_MAX - 1) + 1; ++j)
                  for (k = 0; k < NWAVES; ++k)
                       spectra[iZ][j][k] = spectra[3][j][k] + (spectra[7][j][k] - spectra[3][j][k])/4*(iZ - 3);
        else if (iZ > 7 && iZ < 19)
             for (j = 0; j < NSMOOTH*(SNAP_MAX - 1) + 1; ++j)
                  for (k = 0; k < NWAVES; ++k)
                       spectra[iZ][j][k] = spectra[7][j][k] + (spectra[19][j][k] - spectra[7][j][k])/12*(iZ - 7);
        else if (iZ > 19 && iZ < 39)
             for (j = 0; j < NSMOOTH*(SNAP_MAX - 1) + 1; ++j)
                  for (k = 0; k < NWAVES; ++k)
                       spectra[iZ][j][k] = spectra[19][j][k] + (spectra[39][j][k] - spectra[19][j][k])/20*(iZ - 19);
        else
            ; 
     }
     //


    #if APPARENT || ABSOLUTE
        fp = openfile("mags.dat", "w");
        dl = d_lumin(z0);
        fprintf(fp, "AB magnitudes for galaxies at z = %8.4f, luminosity distance: %10.3E cm\n", z0, dl);
    
        #if APPARENT
            for (k = 0; k < NWAVES; ++k)
                Waveobs[k] = Wave[k]*(1 + z0);

            puts("================ Read the filters ================");
            for(i = 0; i < sizeof(bands)/sizeof(int); ++i)
                if (bands[i] < 8)
                    Thrutot[i] = read_thru(fltlist[bands[i]][1], Thru[i]);
                else if (bands[i] == 8)
                    Thrutot[i] = read_thru36(fltlist[bands[i]][1], Thru[i]);
                else {
                    puts("Error for bands");
                    exit(1);
                }

        #endif
        #if ABSOLUTE
        for(i = 0; i < sizeof(rest_bands)/sizeof(rest_bands[0]); ++i)
            Thrutot_rest[i] = creat_thru(Thru_rest[i], rest_bands[i][0], rest_bands[i][1]);
        #endif

        #if DUSTEXTINC
        //    dust_curve(Wave);
        #endif

        #if LYALPHA
            absorb_table(z0, Waveobs);
        #endif

        #if APPARENT
            for(i = 0; i < sizeof(bands)/sizeof(int); ++i)
                fprintf(fp, "  %5s ", fltlist[bands[i]][0]);
        #endif
        # if ABSOLUTE
        for (i = 0; i < sizeof(rest_title)/sizeof(rest_title[0]); i++)
            fprintf(fp, "  %5s ", rest_title[i]);
        #endif
        fprintf(fp, "\n");

        puts("================== Process the SED of galaxies =============");

        for(n = 0; ;n++) {
            fread(&Nprog, 4, 1, fp0);
            if (Nprog == -1) break;
            if (!(n % 1000))
                printf("Info: Processing galaxy: No. %dk - %dk ...\n", n/1000, n/1000 + 1);

            for (k = 0; k < NWAVES; ++k)
                Energy[k] = 0.0;

            for (i = 0; i < Nprog; i++) {
                fread(&isnap, 4, 1, fp0);
                fread(&metal, 4, 1, fp0);
                fread(&newstars, 4, 1, fp0);

                if (newstars == 0.0) continue;

                if (metal < 1e-3)
                    imetal = 0;
                else if (metal > 4e-2)
                    imetal = 39;
                else
                    imetal= (int)(metal*1000 - 0.5);

                tmp = newstars/NSMOOTH;
                spectrum = spectra[imetal];  //// new added

                if (isnap > 0) {
                for (j = 0; j < NSMOOTH; ++j) {
                    l = isnap*NSMOOTH - j; 
                    for (k = 0; k < NWAVES; ++k)
                        Energy[k] += spectrum[l][k]*tmp;
                }
                }
            
                else
                    for (k = 0; k < NWAVES; ++k) {
                        Energy[k] += spectrum[isnap][k]*newstars;
                    }

            }
        /* Apparent magnitude calculation */
        #if APPARENT
            #if DUSTEXTINC
            #endif
            for (k = 0; k < NWAVES; ++k) {
                Flux[k] = Energy[k]/(4*3.14159*dl*dl*(1 + z0));
            #if DUSTEXTINC
            #endif
            #if LYALPHA
                Flux[k] *= Absorb[k];
            #endif
            }
            for(i = 0; i < sizeof(bands)/sizeof(int); ++i) {
                if (bands[i] < 8)
                    thru_filter(Thru[i], Waveobs, Flux, Flxobs);
                else if (bands[i] == 8)
                    thru_filter36(Thru[i], Waveobs, Flux, Flxobs);
                else {
                    puts("Wrong bands");
                    exit(1);
                }
                fprintf(fp, "%8.3f", abmag(Waveobs, Flxobs, Thrutot[i]));
            }
        #endif
        /* Absolute magnitude calculation */
        #if ABSOLUTE
            for (k = 0; k < NWAVES; ++k)
                Flux[k] = Energy[k]/(4*3.14159*3.08568E19*3.08568E19);
            for (i = 0; i < sizeof(rest_bands)/sizeof(rest_bands[0]); i++) {
                thru_filter02(Thru_rest[i], Wave, Flux, Flxobs);
                fprintf(fp, "%8.3f", abmag(Wave, Flxobs, Thrutot_rest[i]));
            }
        #endif
        fprintf(fp, "\n");
        }
        fclose(fp);
    #else
        puts("Info: No Need for galaxy magnitude calculation?");
    #endif
    fclose(fp0);

    tmend = time(NULL);
    printf("\n Run Finished: %s\n", asctime(localtime(&tmend)));
    printf("Total Run Time: %.0f s\n", difftime(tmend, tmstart));
    puts("==================== END ======================\n");
    return 0;
}


FILE *openfile(char *fname, char *mode)
/* Function: Open the file with specific mode */
{
    FILE *fp;
    if ((fp = fopen(fname, mode)) == NULL) {
        printf("File open error: \"%s\"!\n", fname);
        exit(0);
    }
    printf("File opened: \"%s\"!\n", fname);
    return fp;
}


float read_data(FILE *fp, float *data)
/* Function: read sed from file
     Return: time */
{
    float t;
    int i;
    if (feof(fp)) {
        printf("Error: End of File.\n");
        exit(0);
    }

    fscanf(fp, "%f", &t);
    fseek(fp, -11L, 1);
    for (i = 0; i < NWAVES; ++i) {
        fseek(fp, 30L, 1);
        fscanf(fp, "%f", data++);
        while (fgetc(fp) != '\n');
    }
    return t;
}


/*void dust_curve(float *wave)
// k(lambda) curve 
{
    float lambda, slope1, slope2, intcept1, intcept2;
    int k;
}*/


double d_lumin(float z)
/* Function: Calculate luminosity distance at red shift z. */
{
    double dm, dl;
    float zi, dz;
    dm = 0.0;
    dz = z/10000.0;
    if (Omega_k == 0) {
        for (zi = 0.0; zi <= z; zi += dz)
            dm += dz/sqrt(Omega_M*pow(1 + zi, 3.0) + Omega_L);
        dm *= 2997.9/h0*3.0857E24;
        dl = dm*(1 + z);
    }
    else {
        printf("Error\n");
        exit(0);
    }
    return dl;
}


void absorb_table(float z, float *waveobs)
{
    float tau, zi;
    int k;

    for (k = 0; k < NWAVES; ++k) {
       if (waveobs[k] > 1216.0 && waveobs[k] < 1216.0*(1 + z)) {
           zi = waveobs[k]/1216.0 - 1.0;
           if (zi < 5.5)
               tau = 0.85*pow((1 + zi)/5.0, 4.3);
           else
               tau = 0.15*pow((1 + zi)/5.0, 10.9);
           Absorb[k] = exp(-tau);
       }
       else
           Absorb[k] = 1.0;
    }
}


double read_thru(char *fltname, float *thru)
/* Read filter throughputs */
{
    FILE *openfile(char *, char *);
    FILE *fp;
    double thrutot = 0.0;
    int k;

    fp = openfile(fltname, "r");
    for (k = 0; k < 5; ++k)
        while (fgetc(fp) != '\n');
 
    for (k = 0; k < FILTERUP - FILTERLOW; ++k)
        fscanf(fp, "%*f%*f%f", thru + k);
    close(fp);

    for (k = 0; k < FILTERUP - FILTERLOW; ++k)
         // thrutot += thru[k]*(2.9979E18/(k + FILTERLOW - 1) - 2.9979E18/(k + FILTERLOW));
         thrutot += thru[k]*2.9979E18/((k + FILTERLOW)*(k + FILTERLOW));
    return thrutot;
}


double read_thru36(char *fltname, float *thru)
/* Read filter throughputs */
{
    FILE *openfile(char *, char *);
    FILE *fp;
    double thrutot = 0.0;
    float ktmp, thrutmp;
    int k, i, kint;


    fp = openfile(fltname, "r");
    for (k = 0; k < 2; ++k)
        while (fgetc(fp) != '\n');
    thrutmp = 0.000490;
    k = 0;
    thru[k] = thrutmp;
    for (i = 0; i < 390; ++i) {
        fscanf(fp, "%f%f", &ktmp, &thrutmp);
        kint = (int)ktmp;
        for (; k <= kint - 30810; ++k) {
            thru[k] = thrutmp;
            thrutot += thrutmp*2.9979E18/((k + 30810)*(k + 30810));
        }
    }
    close(fp);

    return thrutot;
}


double creat_thru(float *thru, float bcenter, float bwidth)
/* Create a filter with specific band center and width. */
{
    int k;
    double thrutot = 0.0; 
    for (k = 0; k < FILTERUP - FILTERLOW; ++k)
        if (k >= bcenter - 0.5*bwidth && k <= bcenter + 0.5*bwidth) {
            thru[k] = 1.0;
            thrutot += 2.9979E18/(k*k);
        }
        else
            thru[k] = 0.0;
    return thrutot;
}


void thru_filter(float *thru, float *wave, double *flx, double *flxobs)
/* Spectrum through filters */
{
    int k, n;

    for (k = 0; k < NWAVES; ++k) {
        if ((int)wave[k] >= FILTERLOW && (int)wave[k] < FILTERUP)
            flxobs[k] = flx[k]*thru[(int)wave[k] - FILTERLOW];
        else
            flxobs[k] = 0.0;
    }
}

void thru_filter36(float *thru, float *wave, double *flx, double *flxobs)
/* Spectrum through filters */
{
    int k, n;

    for (k = 0; k < NWAVES; ++k) {
        if ((int)wave[k] >= 30810 && (int)wave[k] < 40103)
            flxobs[k] = flx[k]*thru[(int)wave[k] - 30810];
        else
            flxobs[k] = 0.0;
    }
}


void thru_filter02(float *thru, float *wave, double *flx, double *flxobs)

/* Spectrum through filters */
{
    int k, n;

    for (k = 0; k < NWAVES; ++k) {
        if ((int)wave[k] < FILTERUP - FILTERLOW)
            flxobs[k] = flx[k]*thru[(int)wave[k]];
        else
            flxobs[k] = 0.0;
    }
}


float abmag(float *wave, double *flux, double thrutot)
/* Function: calculate AB magnitude */
{
    double e = 0.0;
    int k;

    for (k = 1; k < NWAVES - 1; ++k)
        e += flux[k]*0.5*(wave[k+1] - wave[k-1]);

    return -2.5*log10(e/thrutot) - 48.6;
}
