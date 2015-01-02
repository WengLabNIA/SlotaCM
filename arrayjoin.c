/* Reading the Agilent gene chip and dye normalization */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>

#define MAXCOL  3000
#define MISSING -9999
#define ONEDYE 1
#define TWODYE 2

typedef struct data_st{
    int ncol;
    int nrow;
    int maxrepl;
    int *n[MAXCOL];
    int *flag[MAXCOL];
    int nSaturated[MAXCOL];
    float *data[MAXCOL];
    float *error[MAXCOL];
    char *header_col[MAXCOL];
    char **featureID;
}DATA;

typedef struct param_st{
    int nrow;
    int nfiles;
    char *dataFiles[MAXCOL];
    char inputFile[300];
    char outputFile[300];
    char aliasFile[300];
    float cv_thresh;
    float error_ratio;
    int dye_swap[MAXCOL];
    int dye_swap_temp[MAXCOL]; // if there exists i: dye_swap[i]==1, dye_swap[i+1]==2
    int array_type;
    int skip_lines;
    int header_line;
    int col_mean[2];
    int col_err[2];
    int col_backgr[2];
    int col_ID;
    int col_condition;
    float condition_min;
    float condition_max;
    char **alias_name;
    char **oligo_name;
    int *oligo;
    int nalias;
    int noligos;
    char *tissue[MAXCOL];
    int replication[MAXCOL];
    int tissueOption;
    int replicOption;
    float maxValue;
    float cutoff;
}PARAM;

void print_line (FILE *fp);
char  *copy_string (char *x);

/***********************************************/
void check (void *x)
/***********************************************/
{
    if (!x){
        printf("Out of memory\n");
        printf("\nHit any key to close the window\n");
        while(!kbhit()){}
        exit(0);
    }
}

/***********************************************/
void error_message (char *message)
/***********************************************/
{
    printf("ERROR: %s\n", message);
    printf("\nHit any key to close the window\n");
    while(!kbhit()){}
    exit(0);
}

/***********************************************/
int  read_to_the_end_of_line (FILE *fp)
/***********************************************/
{
    char a;
    while((a=fgetc(fp))!=EOF && a != '\n');
    if(a==EOF) return(0);
    return(1);
}

/***********************************************/
int count_fields(char *buffer)
/***********************************************/
{
    int i, n=1, len;

    len= strlen(buffer);
    for(i=0; i<len && buffer[i] != '\n'; ++i){
        if(buffer[i] == '\t' && buffer[i+1] !='\n')
            ++n;
    }
    return(n);
}

/***********************************************/
void    split_string (char *line, char **items, int n, int maxlen)
/***********************************************/
{
    char a;
    int j, jpos, pos, len;

    for(j=0; j<n; ++j){
        items[j][0] = '\0';
    }
    len = strlen(line);
    j=0;
    jpos=0;
    pos=0;
    a = '\t';
    while(a != EOF && a != '\n'){
        while(pos<len && (a=line[pos++])!=EOF && a != '\n' && a != '\t'){
            if(jpos<maxlen-1) items[j][jpos++] = a;
        }
        items[j][jpos] = '\0';
        j++;
        jpos=0;
    }
}

/*************************************************************************/
void sortem   (int ie, char **a, int *b)
/*************************************************************************/
{
    int i, j, k, m, p, q, iring;
    int lt[64], ut[64];
    char *ta, *xa;
    int i1, tb, xb;

    j = ie-1;
    m = 1;
    i = 0;

/* If this segment has more than two elements  we split it */
L10:    if ((i1 = j - i - 1) < 0) goto L100;
    else if (i1 == 0) goto L90;
    else goto L15;

/* p is the position of an arbitrary element in the segment we choose the
* middle element. Under certain circumstances it may be advantageous
* to choose p at random. */

L15:
    p = (j + i) / 2;
    ta = a[p];
    a[p] = a[i];
    if(b){
        tb = b[p];
        b[p] = b[i];
    }

L21: /* Start at the beginning of the segment, search for k such that a(k)>t */
    q = j;
    k = i;
L20:    ++k;
    if (k > q) goto L60;
    if (strcmp(a[k],ta)<=0) goto L20;

/* Such an element has now been found now search for a q such that a(q)<t
* starting at the end of the segment. */
L30:    if (strcmp(a[q],ta)<0) goto L40;
    --q;
    if (q > k) goto L30;
    goto L50;

/* a(q) has now been found. we interchange a(q) and a(k) */

L40:    xa = a[k];
    a[k] = a[q];
    a[q] = xa;
    if(b){
        xb = b[k];
        b[k] = b[q];
        b[q] = xb;
    }
/* Update q and search for another pair to interchange: */
    --q;
    goto L20;
L50:    q = k - 1;
L60:
/* The upwards search has now met the downwards search: */
    a[i] = a[q];
    a[q] = ta;
    if(b){
        b[i] = b[q];
        b[q] = tb;
    }
L65:
/* The segment is now divided in three parts: (i,q-1),(q),(q+1,j) */
/* store the position of the largest segment in lt and ut */
    if (q << 1 <= i + j) goto L70;
    lt[m - 1] = i;
    ut[m - 1] = q - 1;
    i = q + 1;
    goto L80;
L70:    lt[m - 1] = q + 1;
    ut[m - 1] = j;
    j = q - 1;
/* Update m and split the new smaller segment */
L80:    ++m;
    goto L10;

/* We arrive here if the segment has  two elements we test to see if */
/* the segment is properly ordered if not, we perform an interchange */
L90:
    if (strcmp(a[i],a[j])<=0) goto L100;
    xa = a[i];
    a[i] = a[j];
    a[j] = xa;
    if(b){
        xb = b[i];
        b[i] = b[j];
        b[j] = xb;
    }
L95:

/* If lt and ut contain more segments to be sorted repeat process: */
L100:    --m;
     if (m <= 0) goto L110;
     i = lt[m - 1];
     j = ut[m - 1];
     goto L10;
L110:    return;
} /* sortem_ */

/***********************************************/
int  find_string (char *word, char **a, int n)
/***********************************************/
{
    int i, j=-1, sign;

    if(n<=0) return(-1);
    i = n/2;
    while((sign=strcmp(word,a[i]))){
        if(sign>0){
            if(n-i<=1){ return(-1); }
            j = i;
            i = (n+i)/2;
        }
        else{
            if(i-j<=1){ return(-1); }
            n = i;
            i = (j+n)/2;
        }
    }
    return(i);
}

/***********************************************/
DATA *read_data (PARAM *p)
/***********************************************/
{
    FILE *inputFile;
    DATA *d;
    int i, j, ifile, irow, icol, irow1, ncol_data, idye, newrow, flag, len;
    int iline;
    float mean[2], error[2], backgr[2], swap, value;
    char *line, *filename, oligo_name[40], data_name[200], buffer[50], *s;
    char **items;
    int *index, found, last_sort=0;
    char **oligo_list;
    int *oligo_index;

    line = (char*)malloc(4000*sizeof(char));
    d = (DATA*)calloc(1,sizeof(DATA));
    d->ncol = p->nfiles*p->array_type;

    inputFile = fopen(p->dataFiles[0], "r");
    if(!inputFile){
        printf("ERROR: File %s not found\n", p->dataFiles[0]);
        while(!kbhit());
        exit(0);
    }
    iline = 0;
    if(p->header_line >= 0){
        for(;iline<p->header_line;++iline)
            read_to_the_end_of_line(inputFile);
        fgets(line,3999,inputFile);
        ncol_data = count_fields(line);
        items = (char**)malloc(ncol_data*sizeof(char*));
        for(i=0;i<ncol_data;++i){
            items[i] = (char*)malloc(100*sizeof(char));
        }
        split_string(line, items, ncol_data,100);
        for(i=0;i<ncol_data;++i){
            if(!strcmp(items[i],"ControlType") && p->col_condition<0) p->col_condition = i;
            else if(!strcmp(items[i],"ControlType") && p->col_condition<0) p->col_condition = i;
            else if(!strcmp(items[i],"ProbeName") && p->col_ID<0) p->col_ID = i;
            else if(!strcmp(items[i],"gProcessedSignal") && p->col_mean[0]<0) p->col_mean[0] = i;
            else if(!strcmp(items[i],"rProcessedSignal") && p->col_mean[1]<0) p->col_mean[1] = i;
            else if(!strcmp(items[i],"gProcessedSigError") && p->col_err[0]<0) p->col_err[0] = i;
            else if(!strcmp(items[i],"rProcessedSigError") && p->col_err[1]<0) p->col_err[1] = i;
        }
    }
    if(p->col_ID < 0 || p->col_mean[0] < 0 || p->col_mean[1] <0){
        error_message("Column headers not found");
    }
    printf("Column number for feature ID: %d\n",p->col_ID+1);
    printf("Column number for 1st dye mean: %d\n",p->col_mean[0]+1);
    if(p->array_type==TWODYE){
        printf("Column number for 2nd dye mean: %d\n",p->col_mean[1]+1);
    }
    if(p->col_err[0]>=0) printf("Column number for 1st dye error: %d\n",p->col_err[0]+1);
    if(p->array_type==TWODYE){
        if(p->col_err[1]>=0) printf("Column number for 2nd dye error: %d\n",p->col_err[1]+1);
    }
    if(p->col_backgr[0]>=0) printf("Column# for 1st dye background: %d\n",p->col_backgr[0]+1);
    if(p->array_type==TWODYE){
        if(p->col_backgr[1]>=0) printf("Column# for 2nd dye background: %d\n",p->col_backgr[1]+1);
    }
    if(p->col_condition>=0){
        printf("Column number for condition that row is included: %d\n",p->col_condition+1);
        printf("Minimum condition value: %f\n", p->condition_min);
        printf("Maximum condition value: %f\n", p->condition_max);
    }
    for(;iline<p->skip_lines;++iline)
        read_to_the_end_of_line(inputFile);
    while(fgets(line,3999,inputFile) && strlen(line)>=2){
        ++p->nrow;
    }
    fclose(inputFile);
    d->nrow = p->nrow;
    if(p->col_mean[0] > ncol_data)
        error_message("Column number for mean1 is too large");
    if(p->array_type==TWODYE && p->col_mean[1] > ncol_data)
        error_message("Column number for mean2 is too large");
    if(p->col_err[0]>=0 && p->col_err[0] >= ncol_data)
        error_message("Column number for error1 is too large");
    if(p->array_type==TWODYE && p->col_err[1]>=0 && p->col_err[1] >= ncol_data)
        error_message("Column number for mean2 is too large");
    if(p->col_ID >= ncol_data)
        error_message("Column number for gene ID is too large");
    if(p->col_condition>=0 && p->col_condition >= ncol_data)
        error_message("Column number for condition is too large");

    oligo_list = (char**)malloc(d->nrow*sizeof(char*));
    oligo_index = (int*)malloc(d->nrow*sizeof(int));
    index = (int*)malloc(d->nrow*sizeof(int));
    d->featureID = (char**)malloc(d->nrow*sizeof(char*));
    for(irow=0;irow<d->nrow;++irow){
        oligo_index[irow] = irow;
    }
    for(i=0;i<d->ncol;++i){
        d->n[i] = (int*)calloc(d->nrow,sizeof(int));
        d->flag[i] = (int*)calloc(d->nrow,sizeof(int));
        d->data[i] = (float*)malloc(d->nrow*sizeof(float));
        d->error[i] = (float*)malloc(d->nrow*sizeof(float));
        for(irow=0;irow<d->nrow;++irow){
            d->data[i][irow] = MISSING;
            d->error[i][irow] = MISSING;
        }
    }
    for(ifile=0; ifile<p->nfiles; ++ifile){
        s = p->dataFiles[ifile];
        len = strlen(s);
        filename = strrchr(s,'\\');
        if(filename < s || filename > s+len) filename = strrchr(s,'/');
        if(filename < s || filename > s+len) filename = s;
        else filename++;
        inputFile = fopen(s, "r");
        if(!inputFile){
            printf("Input file %s not found\nHit any key to close the window\n", filename);
            while(!kbhit()){}
            exit(0);
        }
        printf("Reading file %s\n", filename);
        if(p->tissueOption){
            strcpy(data_name,p->tissue[ifile]);
            if(p->replicOption){
                strcat(data_name,"_rep");
                sprintf(buffer,"%d",p->replication[ifile]);
                strcat(data_name,buffer);
            }
        }else{
            strcpy(data_name, filename);
        }
        for(idye=0; idye<p->array_type; ++idye){
            icol = ifile*p->array_type+idye;
            if(idye) strcpy(data_name,"Reference");
            d->header_col[icol] = copy_string(data_name);
        }
        for(i=0;i<p->skip_lines;++i){
            read_to_the_end_of_line(inputFile);
        }
        if(ifile == 0) newrow = 0;
        for(irow=0;irow<d->nrow;++irow){
            fgets(line,2000,inputFile);
            if(ncol_data != count_fields(line)){
                printf("Error: Wrong number of fields in row %d\n", irow);
                printf("\nHit any key to close the window\n");
                while(!kbhit());
                exit(0);
            }
            split_string(line, items, ncol_data,70);
            if(p->col_condition >= 1){
                sscanf(items[p->col_condition],"%f", &value);
                if(value < p->condition_min || value > p->condition_max)
                    continue;
            }
            if(ifile == 0){
                if(!strlen(items[p->col_ID])) error_message("Empty feature ID field");
                /* Find oligo name from alias */

                found = find_string(items[p->col_ID],p->alias_name,p->nalias);
                if(found >= 0) strcpy(oligo_name,p->oligo_name[p->oligo[found]]);
                else strcpy(oligo_name,items[p->col_ID]);

                found = find_string(oligo_name,oligo_list,last_sort);
                for(i=last_sort; i<newrow && found==-1; ++i){
                    if(!strcmp(oligo_name,oligo_list[i])) found=i;
                }
                if(found >=0) index[irow] = oligo_index[found];
                else{
                    index[irow] = newrow;
                    d->featureID[newrow] = copy_string(oligo_name);
                    oligo_list[newrow] = d->featureID[newrow];
                    ++newrow;
                    if(newrow%1000==0){
                        sortem(newrow,oligo_list,oligo_index);
                        last_sort = newrow;
                    }
                }
            }
            for(idye=0; idye<p->array_type; ++idye){
                if(sscanf(items[p->col_mean[idye]],"%f",&mean[idye])!=1)error_message("Cannot read mean value");
                if(p->col_err[idye] >= 0){
                    if(sscanf(items[p->col_err[idye]],"%f",&error[idye])!=1)error_message("Cannot read error value");
                }else{
                    error[idye] = 0.00001;
                }
                if(p->col_backgr[idye] >= 0){
                    if(sscanf(items[p->col_backgr[idye]],"%f",&backgr[idye])!=1)error_message("Cannot read background value");
                }else{
                    backgr[idye] = 0;
                }
                mean[idye] -= backgr[idye];
            }
            if(p->dye_swap[ifile] == 1){
                swap = mean[0];
                mean[0] = mean[1];
                mean[1] = swap;
                swap = error[0];
                error[0] = error[1];
                error[1] = swap;
            }
            irow1 = index[irow];
            for(idye=0; idye<p->array_type; ++idye){
                icol = ifile*p->array_type+idye;
                flag = 0;
                if(mean[idye] <= 0 || error[idye]/mean[idye] >= p->cv_thresh)
                    flag = 1;

                if(d->n[icol][irow1]==0){
                    d->data[icol][irow1] = mean[idye];
                    d->error[icol][irow1] = error[idye];
                    if(flag)
                        d->flag[icol][irow1] = 1;
                    d->n[icol][irow1]++;
                }
                else{
                    if(d->flag[icol][irow1] && flag || !d->flag[icol][irow1] && !flag){
                        d->data[icol][irow1] += mean[idye];
                        d->error[icol][irow1] += error[idye];
                        d->n[icol][irow1]++;
                    }
                    else if(d->flag[icol][irow1]){
                        d->data[icol][irow1] = mean[idye];
                        d->error[icol][irow1] = error[idye];
                        d->flag[icol][irow1] = 0;
                        d->n[icol][irow1] = 1;
                    }
                    else continue;
                }
            }
        }
        fclose(inputFile);
        for(idye=0; idye<p->array_type; ++idye){
            icol = ifile*p->array_type+idye;
            for(irow=0;irow<d->nrow;++irow){
                if(d->n[icol][irow]>1){
                    d->data[icol][irow] /= d->n[icol][irow];
                    d->error[icol][irow] /= d->n[icol][irow];
                }
                if(d->data[icol][irow] > p->maxValue){
printf("%d\t%d\t%f\n",icol,irow,d->data[icol][irow]);
                    d->nSaturated[icol]++;
                }
            }
            printf("%d saturated values in column %d\n", d->nSaturated[icol], icol+1);
        }
    }
    d->nrow = newrow;
    p->nrow = d->nrow;
    printf ("Rows = %d\n",p->nrow);
    free(line);
    free(index);
    free(oligo_list);
    free(oligo_index);
    if(p->header_line >= 0){
        for(i=0;i<ncol_data;++i)
            free(items[i]);
        free(items);
    }
    return(d);
}

/***********************************************/
void    remove_bad_data  (DATA *d, PARAM *p)
/***********************************************/
{
    float err_func[300];
    int n_err_func[300];
    int icol, irow, i, NERR, nmissing, nsurrogates, nsaturated, saturated;
    float xmax, dx, error;

    nmissing = 0;
    nsurrogates = 0;
    nsaturated = 0;
    NERR = d->nrow/20;
    if(NERR > 300) NERR = 300;
    for(icol=0; icol<d->ncol; icol++){

        /* Estimate the error function */
        xmax = 0;
        for(irow=0; irow<d->nrow; ++irow){
            if(d->data[icol][irow] > MISSING && xmax < d->data[icol][irow])
                xmax = d->data[icol][irow];
        }
        dx = sqrt(xmax)/NERR;
        for(i=0; i<NERR; ++i){
            err_func[i] = 0;
            n_err_func[i] = 0;
        }
        for(irow=0; irow<d->nrow; ++irow){
            if(d->error[icol][irow] <= MISSING || d->error[icol][irow]/d->data[icol][irow]>p->cv_thresh)
                continue;
            if(d->data[icol][irow] <= 1.0e-10)
                i = 0;
            else
                i = sqrt(d->data[icol][irow])/dx;
            if(i>=NERR) i = NERR-1;
            err_func[i] += d->error[icol][irow];
            n_err_func[i]++;
        }
        for(i=0; i<NERR; ++i){
            if(n_err_func[i] > 3 || (i==0 && n_err_func[i] > 0))
                err_func[i] /= n_err_func[i];
            else if(i>0 && n_err_func[i-1] > 3){
                err_func[i] = err_func[i-1];
                n_err_func[i] = n_err_func[i-1];
            }
            else err_func[i] = -1;
        }
        for(i=NERR-1; i>=0; --i){
            if(err_func[i]==-1 && i<NERR-1)
                err_func[i] = err_func[i+1];
        }
        saturated = 0;
        if(d->nSaturated[icol] > 0.01*d->nrow){ saturated = 1; }

        /* Remove/fix bad data */
        for(irow=0; irow<d->nrow; ++irow){
            if(saturated && d->data[icol][irow]>p->maxValue){
                d->data[icol][irow] = MISSING;
                ++nsaturated;
                continue;
            }
            if(p->cv_thresh <= 0 || !d->flag[icol][irow])
                continue;
            if(d->data[icol][irow] <= 1.0e-10)
                i = 0;
            else
                i = sqrt(d->data[icol][irow])/dx;
            error = err_func[i];
            if(d->error[icol][irow]/error > p->error_ratio){
                d->data[icol][irow] = MISSING;
                ++nmissing;
            }
            else if (d->data[icol][irow] < error){
                d->data[icol][irow] = error;
                ++nsurrogates;
            }
        }
    }
    if(nmissing)
        printf("Data with relative error > %f: %d values replaced by -9999\n", p->cv_thresh, nmissing);
    if(nsaturated)
        printf("Data with saturated values > %f: %d values replaced by -9999\n", p->maxValue, nsaturated);
    if(nsurrogates)
        printf("%d low data numbers were substituted with surrogates (=error value)\n", nsurrogates);
}

/***********************************************/
void print_output  (char *filename, DATA *d)
/***********************************************/
{
    FILE *outputFile;
    int icol,irow;

    outputFile = fopen(filename, "w");
    if(!outputFile){
        printf("Output file %s not writable!\n", filename);
        printf("\nHit any key to close the window\n");
        while(!kbhit()){}
        exit(0);
    }

    fprintf(outputFile,"FeatureID\t");
    for(icol=0;icol<d->ncol-1;++icol)
        fprintf(outputFile,"%s\t", d->header_col[icol]);
    fprintf(outputFile,"%s\n", d->header_col[d->ncol-1]);

    for(irow=0;irow<d->nrow;++irow){
        fprintf(outputFile,"%s\t",d->featureID[irow]);
        for(icol=0;icol<d->ncol-1;++icol){
            fprintf(outputFile,"%.4f\t", d->data[icol][irow]);
        }
        fprintf(outputFile,"%.4f\n", d->data[d->ncol-1][irow]);
    }
    fclose(outputFile);
}

/**********************************************/
float      read_float   (FILE *fp, char *description)
/**********************************************/
{
float x;

    if (fscanf (fp,"%f",&x) != 1)
    {
        print_line(fp);
        printf("Error reading parameter: %s\n", description);
        while(!kbhit());
        exit(0);
    }
    printf("%f\t%s\n", x, description);
    return (x);
}

/****************************/
void print_line (FILE *fp)
/***************************/
{
    char buffer[21];
    fgets(buffer, 20, fp);
    buffer[20] = '\0';
    printf("%s\n", buffer);
}

/**********************************************/
int   read_int  (FILE *fp, char *description)
/**********************************************/
{
    int x, index;

    if (fscanf (fp,"%d",&x) != 1)
    {
        print_line(fp);
        printf("Error reading parameter: %s\n", description);
        while(!kbhit());
        exit(0);
    }
    index = fgetc (fp);
    if (index == '.')   /* Float number is read instead of integer */
    {
        print_line(fp);
        printf("Integer is expected, in parameter %s", description);
        while(!kbhit());
        exit(0);
    }
    else
        ungetc (index, fp);
    printf("%d\t%s\n", x, description);
    return (x);
}

/***********************************************/
char  *copy_string (char *x)
/***********************************************/
{
    int len;
    char *y;

    len = strlen(x);
    check(y = (char*)malloc((len+1)*sizeof(char)));
    strcpy(y, x);
    return(y);
}

/***********************************************/
PARAM *read_parameters (char **argv, int nargs)
/***********************************************/
{
    FILE *fp;
    PARAM *p;
    int i, j, len;
    char buffer[300], **items;
    int iarg=1, ialias, ioligo, size[4]={3,200,200,15};


    p = (PARAM*)calloc(1,sizeof(PARAM));
    p->cv_thresh = 1000000000;
    p->error_ratio = 3;
    p->array_type = TWODYE;
    p->skip_lines = 10;
    p->header_line = 10;
    p->maxValue = 1.0E20;
    p->cutoff = MISSING;
    while(iarg < nargs){
        if(!strcmpi(argv[iarg],"-i")) strcpy(p->inputFile,argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-o")) strcpy(p->outputFile,argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-a")) strcpy(p->aliasFile,argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-t1")) sscanf(argv[++iarg],"%f",&p->cv_thresh);
        else if(!strcmpi(argv[iarg],"-t2")) sscanf(argv[++iarg],"%f",&p->error_ratio);
        else if(!strcmpi(argv[iarg],"-d")) p->array_type = atoi(argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-r")) p->skip_lines = atoi(argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-h")) p->header_line = atoi(argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-m")){
            if(sscanf(argv[++iarg],"%d,%d",&p->col_mean[0],&p->col_mean[1])==1) p->array_type = ONEDYE;
        }
        else if(!strcmpi(argv[iarg],"-e")) sscanf(argv[++iarg],"%d,%d",&p->col_err[0],&p->col_err[1]);
        else if(!strcmpi(argv[iarg],"-bg")) sscanf(argv[++iarg],"%d,%d",&p->col_backgr[0],&p->col_backgr[1]);
        else if(!strcmpi(argv[iarg],"-id")) p->col_ID = atoi(argv[++iarg]);
        else if(!strcmpi(argv[iarg],"-c")) sscanf(argv[++iarg],"%d,%f,%f",&p->col_condition,&p->condition_min,&p->condition_max);
        else if(!strcmpi(argv[iarg],"-max")) sscanf(argv[++iarg],"%f",&p->maxValue);
        else if(!strcmpi(argv[iarg],"-cut")) sscanf(argv[++iarg],"%f",&p->cutoff);
        else{
            printf("Wrong option %s\n", argv[iarg]);
            while(!kbhit());
            exit(0);
        }
        ++iarg;
    }
    if(!strlen(p->inputFile)){
        printf("Input file name: ");
        scanf("%s", p->inputFile);
    }
    if(!strlen(p->outputFile)){
        printf("Output file name: ");
        scanf("%s", p->outputFile);
    }
    if(p->condition_max < p->condition_min){
        printf("Condition minimun should not be greater than condition maximum\n");
        while(!kbhit());
        exit(0);
    }
    printf("Input file: %s\n", p->inputFile);
    printf("Output file: %s\n", p->outputFile);
    if(p->aliasFile) printf("Alias file: %s\n", p->aliasFile);
    printf("Output file: %s\n", p->outputFile);
    printf("Threshold for acceptable (Error/Mean) ratio: %f\n",p->cv_thresh);
    printf("Threshold for (feature error/expected error): %f\n",p->error_ratio);
    printf("Maximum value (saturation): %f\n",p->maxValue);
    printf("Number of dyes in the array: %d\n", p->array_type);
    printf("Number of header lines to skip: %d\n",p->skip_lines);
    printf("Line number with headers: %d\n",p->header_line);
    if(p->col_mean[0]) printf("Column# for 1st dye mean: %d\n",p->col_mean[0]);
    if(p->array_type==TWODYE){
        if(p->col_mean[1]) printf("Column# for 2nd dye mean: %d\n",p->col_mean[1]);
    }
    if(p->col_err[0]) printf("Column# for 1st dye error: %d\n",p->col_err[0]);
    if(p->array_type==TWODYE){
        if(p->col_err[1]) printf("Column# for 2nd dye error: %d\n",p->col_err[1]);
    }
    if(p->col_backgr[0]) printf("Column# for 1st dye background: %d\n",p->col_backgr[0]);
    if(p->array_type==TWODYE){
        if(p->col_backgr[1]) printf("Column# for 2nd dye background: %d\n",p->col_backgr[1]);
    }

    if(p->col_ID) printf("Column number for feature ID: %d\n",p->col_ID);
    if(p->col_condition){
        printf("Column number for condition that row is included: %d\n",p->col_condition);
        printf("Minimum condition value: %f\n", p->condition_min);
        printf("Maximum condition value: %f\n", p->condition_max);
    }
    --p->header_line;
    --p->col_mean[0];
    --p->col_mean[1];
    --p->col_err[0];
    --p->col_err[1];
    --p->col_backgr[0];
    --p->col_backgr[1];
    --p->col_ID;
    --p->col_condition;

    /* Read input file */
    fp = fopen(p->inputFile,"r");
    if(!fp){ printf("Input file %s not found\n",p->inputFile); exit(0); }
    p->nfiles=0;
    items = (char**)malloc(4*sizeof(char*));
    for(i=0;i<4;++i){
        items[i] = (char*)malloc(size[i]*sizeof(char));
    }
    p->tissueOption = 1;
    p->replicOption = 1;
    while(fgets(buffer,299,fp)){
        if(strlen(buffer)<3) continue;
        split_string(buffer, items, 4, 200);
        i = p->nfiles;
        if(!strlen(items[1])) continue;
        p->dye_swap[i] = atoi(items[0]);
        if(p->dye_swap[i] < 0 || p->dye_swap[i] > 2)
            error_message("ERR: Dye swap indicator should be either 0, 1 or 2\n");
        if(p->dye_swap[0] == 2)
            error_message("First file in a list cannot have type = 2");

        p->dataFiles[i] = copy_string(items[1]);
        if(strlen(items[2])){
            p->tissue[i] = copy_string(items[2]);
        }else{
            p->tissueOption = 0;
            p->replicOption = 0;
        }
        if(strlen(items[3])){
            p->replication[i] = atoi(items[3]);
        }else{
            p->replicOption = 0;
        }
        p->nfiles++;
    }
    fclose(fp);
    if(!p->nfiles){
        printf("Incorrect parameter file. See readme.txt file for the format\n");
        exit(0);
    }

    /* Read alias file */
    if(!strlen(p->aliasFile)){
        for(i=0;i<4;++i) free(items[i]);
        free(items);
        return(p);
    }
    fp = fopen(p->aliasFile,"r");
    if(!fp) error_message("Alias file you entered not found");
    p->nalias=0;
    p->noligos=0;
    while(fgets(buffer,200,fp)){
        len = strlen(buffer);
        if(len > 2){ ++p->nalias; ++p->noligos; }
        for(i=0;i<len;++i){
            if(buffer[i] == ',') ++p->nalias;
        }
    }
    rewind(fp);
    if(!p->nalias){
        printf("Empty alias file\n");
        exit(0);
    }

    p->alias_name = (char**)malloc(p->nalias*sizeof(char*));
    p->oligo_name = (char**)malloc(p->noligos*sizeof(char*));
    p->oligo = (int*)malloc(p->nalias*sizeof(int));
    ialias=0;
    ioligo=0;
    while(fgets(buffer,200,fp)){
        len = strlen(buffer);
        for(i=0;i<i<len && buffer[i] != '\t';++i){}
        if(buffer[i] != '\t') continue;
        buffer[i] = '\0';
        p->oligo_name[ioligo] = copy_string(buffer);
        ++i;
        while(i<len){
            j = i;
            for(;i<len && buffer[i] != ',' && buffer[i] != '\n'&& buffer[i] != '\t'&& buffer[i] != ' ';++i){}
            if(i < len) buffer[i] = '\0';
            if(i > j){
                p->alias_name[ialias] = copy_string(&buffer[j]);
                p->oligo[ialias] = ioligo;
                ++ialias;
            }
            ++i;
        }
        ++ioligo;
    }
    fclose(fp);
    sortem(p->nalias,p->alias_name,p->oligo);
    for(i=0;i<4;++i) free(items[i]);
    free(items);
    return(p);
}

/**********************************************/
float average (DATA *d, int iRow, int iStart, int iFinish, int arr_type, int if_ref)
/**********************************************/
{
    float sum;
    int i, number, current;

    number = 0;
    sum = 0.0;
    for (i = iStart; i < iFinish; i++) {
        current = arr_type * i + if_ref;
        if (d->data[current][iRow] > MISSING) {
            sum += d->data[current][iRow];
            number++;
        }
    }
    if (number > 0) {
        return (sum / number);
    } else {
        return (MISSING);
    }

}

/***********************************************/
DATA *scan_averaging (DATA *d, PARAM *p)
/***********************************************/
{
    DATA *d1;
    int nScanCopies = 0, i, iFile, iCol, iRow, iStart, iFinish, iZero;

    // SPECIFY:   d->ncol, d->nrow, d->header_col, d->featureID, d->data
    for(iFile = 0; iFile < p->nfiles; ++iFile){
        if(p->dye_swap[iFile] == 2)
            nScanCopies++;
    }
    if (!nScanCopies)
        return(d);

    d1 = (DATA*)calloc(1,sizeof(DATA));
    d1->nrow = d->nrow;
    d1->ncol = (p->nfiles - nScanCopies) * p->array_type;
    d1->featureID = (char**)malloc(d1->nrow*sizeof(char*));
    for(i = 0; i < d1->nrow; ++i) {
        d1->featureID[i] = copy_string(d->featureID[i]);
    }
    for(i = 0; i < d1->ncol; ++i) {
        d1->data[i] = (float*)malloc(d1->nrow*sizeof(float));
    }
    iCol = 0;
    for (iFile = 0; iFile < p->nfiles; ++iFile) {
        if (p->dye_swap[iFile] != 2) {
            d1->header_col[iCol++] = copy_string(d->header_col[iFile * p->array_type]);
            if(p->array_type == TWODYE) {
                d1->header_col[iCol++] = copy_string("Reference");
            }
        }
    }
    for(iRow = 0; iRow < d1->nrow; ++iRow) {
        iStart = 0;
        for (iZero = 0; iZero < (p->nfiles - nScanCopies); ++iZero) {
            iCol = iZero * p->array_type;
            iFinish = iStart + 1;
            while (p->dye_swap[iFinish] == 2)  iFinish++;
            d1->data[iCol][iRow] = average(d, iRow, iStart, iFinish, p->array_type, 0);
            if (p->array_type == TWODYE) {
                d1->data[iCol+1][iRow] = average(d, iRow, iStart, iFinish, p->array_type, 1);
            }
            iStart = iFinish;
        }
    }
    return(d1);
}

/***********************************************/
void dye_swap_processing_before (PARAM *p)
/***********************************************/
{
    int iFile, i;

    for (iFile = 0; iFile < p->nfiles; iFile++)
        p->dye_swap_temp[iFile] = p->dye_swap[iFile];

    // every sequence '12...2,0 or 1' change to 11...1,0 or 1
    i = 1;
    while (i < p->nfiles) {
        if (p->dye_swap[i-1]==1 && p->dye_swap[i]==2) {
            while (i < p->nfiles && p->dye_swap[i]==2) {
                p->dye_swap[i] = 1;
                i++;
            }
        }
        i++;
    }
}

/***********************************************/
void dye_swap_processing_after (PARAM *p)
/***********************************************/
{
    int iFile;

    for (iFile = 0; iFile < p->nfiles; iFile++)
        p->dye_swap[iFile] = p->dye_swap_temp[iFile];
}

/***********************************************/
void cut_off (DATA *d, PARAM *p)
/***********************************************/
{
int irow, icol;

if(p->cutoff <= MISSING) return;
if(p->cutoff <= 0){
    for(irow = 0; irow < d->nrow; irow++){
        for(icol = 0; icol < d->ncol; icol++){
            if (d->data[icol][irow] < p->cutoff) d->data[icol][irow] = p->cutoff;
        }
    }
}
else{
    float range = p->cutoff*0.9;
    for(irow = 0; irow < d->nrow; irow++){
        for(icol = 0; icol < d->ncol; icol++){
            float x = d->data[icol][irow];
            if (x < p->cutoff)
                d->data[icol][irow] = p->cutoff + range*(exp((x-p->cutoff)/range)-1)*rand()/RAND_MAX;
        }
    }
}
return;
}

/***********************************************/
int main (int argc, char **argv)
/***********************************************/
{
    DATA *d, *d1;
    PARAM *p;

    printf("Syntax: arrayjoin [-i input][-o output][-a aliasFile][-t1 error/mean][-t2 error/avr.error][-d dye#][-r rowsSkip][-h headerRow][-m meanCols][-e errorCols][-id idCol][-c col,min,max][-cut cutoff]\n");
    p = read_parameters(argv,argc);

    dye_swap_processing_before(p);

    d = read_data(p);
    remove_bad_data(d, p);

    dye_swap_processing_after(p);

    d1 = scan_averaging(d, p);
    cut_off(d1, p);

    print_output(p->outputFile, d1);

    printf("\nHit any key to close the window\n");
    while(!kbhit()){}
    return(0);
}

