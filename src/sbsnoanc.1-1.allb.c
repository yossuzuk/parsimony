#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "sbsnoanc.h"

int main(int argc, char *argv[])
{

	char seq1[MAXSEQ];
	char seq2[MAXSEQ];

	struct codons{char aa; char property;};
	struct codons codon[5][5][5];

/* hold information on daughter node */
/******
        struct pair {char name1[20]; char name2[20]; double blen; int c1[5]; int c2[5]; int c3[5]; char aa[21];};
******/
        struct pair {char name1[20]; char name2[20]; double blen; char c1[5]; char c2[5]; char c3[5]; char aa[21];};
        struct pair pairs[MAXPAIR];

	FILE *fin;
	char buffer[MAXBUF];
/******
	char sequon[MAXSEQ];
	char bufnodea1[MAXDIG], bufnodea2[MAXDIG];
******/
	char oya[6], ko[6];
	int aaa, aab, aac, aad, npairs;
	int baa, bab;
	int sitenum, entnum, position;
	int site1, site2, site3, site4, site5, site6, site7, site8, site9;
	int site11, site12, site13, site14, site15, site16, site17, site18, site19;
	int site21, site22, site23, site24, site25, site26, site27, site28, site29;
	int oyacodon[4];


	double site, diff;
	double seqsyn, seqinv, seqsumlen, nonsyn, noninv, nonsumlen;
	double intseqsyn, intseqinv, intseqsumlen, intnonsyn, intnoninv, intnonsumlen;
	double wsynsite, wnonsite, tsynchge, tnonchge;
	char cod1[4], cod2[4];
	double codres[9];
	double branlen;

	char accept[2000][21];
	char flag;
	char atta1, atta2;
	int zaa, zab, zac, zad;
	double wconsite, tconchge;

	char blena[20];
	int brannum;

/******
	struct kodomos {int c1[5]; int c2[5]; int c3[5]; char aa[21];};
******/

	struct kodomos {char c1[5]; char c2[5]; char c3[5]; char aa[21];};
	struct kodomos kodomo[3];
	int saa, sab, paa, pab, kaa, kab, kac;
	int maa;
	int numc1, numc2, numc3, numaa;
	int konum;
	char seqatta, onajiattac1, onajiattac2, onajiattac3, onajiattaaa, aminoatta, bufamino;
	char bufname[20];
	char bufseq[2000];
	double chgnumc1, chgnumc2, chgnumc3, chgnumaa;


	if(argc!=5)
	{
		printf("sbsnoanc.x : a program of adaptsite-p eliminating structural constraints\n"); fflush(stdout);
		printf("use        : sbsstr.x extseqfile oyakofile acceptableaminoacids kappa\n"); fflush(stdout);
		printf("example    : sbsstr.x extseqfile oyakofile acceptable 6.0\n"); fflush(stdout);
		exit(1);
	}



/******
	for(aaa=0; aaa<MAXSEQ; aaa++)
	{
		sequon[aaa]='N';
	}
******/



/* mutation rate */
	mu[1][2]=atof(argv[4])*1.0;
	mu[1][3]=1.0;
	mu[1][4]=1.0;
	mu[2][1]=atof(argv[4])*1.0;
	mu[2][3]=1.0;
	mu[2][4]=1.0;
	mu[3][1]=1.0;
	mu[3][2]=1.0;
	mu[3][4]=atof(argv[4])*1.0;
	mu[4][1]=1.0;
	mu[4][2]=1.0;
	mu[4][3]=atof(argv[4])*1.0;

	codon[1][1][1].aa='F'; codon[1][1][1].property='E'; codon[1][1][2].aa='F'; codon[1][1][2].property='E';
	codon[1][1][3].aa='L'; codon[1][1][3].property='E'; codon[1][1][4].aa='L'; codon[1][1][4].property='E';
	codon[2][1][1].aa='L'; codon[2][1][1].property='E'; codon[2][1][2].aa='L'; codon[2][1][2].property='E';
	codon[2][1][3].aa='L'; codon[2][1][3].property='E'; codon[2][1][4].aa='L'; codon[2][1][4].property='E';
	codon[3][1][1].aa='I'; codon[3][1][1].property='E'; codon[3][1][2].aa='I'; codon[3][1][2].property='E';
	codon[3][1][3].aa='I'; codon[3][1][3].property='E'; codon[3][1][4].aa='M'; codon[3][1][4].property='E';
	codon[4][1][1].aa='V'; codon[4][1][1].property='E'; codon[4][1][2].aa='V'; codon[4][1][2].property='E';
	codon[4][1][3].aa='V'; codon[4][1][3].property='E'; codon[4][1][4].aa='V'; codon[4][1][4].property='E';
	codon[1][2][1].aa='S'; codon[1][2][1].property='E'; codon[1][2][2].aa='S'; codon[1][2][2].property='E';
	codon[1][2][3].aa='S'; codon[1][2][3].property='E'; codon[1][2][4].aa='S'; codon[1][2][4].property='E';
	codon[2][2][1].aa='P'; codon[2][2][1].property='E'; codon[2][2][2].aa='P'; codon[2][2][2].property='E';
	codon[2][2][3].aa='P'; codon[2][2][3].property='E'; codon[2][2][4].aa='P'; codon[2][2][4].property='E';
	codon[3][2][1].aa='T'; codon[3][2][1].property='E'; codon[3][2][2].aa='T'; codon[3][2][2].property='E';
	codon[3][2][3].aa='T'; codon[3][2][3].property='E'; codon[3][2][4].aa='T'; codon[3][2][4].property='E';
	codon[4][2][1].aa='A'; codon[4][2][1].property='E'; codon[4][2][2].aa='A'; codon[4][2][2].property='E';
	codon[4][2][3].aa='A'; codon[4][2][3].property='E'; codon[4][2][4].aa='A'; codon[4][2][4].property='E';
	codon[1][3][1].aa='Y'; codon[1][3][1].property='E'; codon[1][3][2].aa='Y'; codon[1][3][2].property='E';
	codon[1][3][3].aa='X'; codon[1][3][3].property='X'; codon[1][3][4].aa='X'; codon[1][3][4].property='X';
	codon[2][3][1].aa='H'; codon[2][3][1].property='P'; codon[2][3][2].aa='H'; codon[2][3][2].property='P';
	codon[2][3][3].aa='Q'; codon[2][3][3].property='E'; codon[2][3][4].aa='Q'; codon[2][3][4].property='E';
	codon[3][3][1].aa='N'; codon[3][3][1].property='E'; codon[3][3][2].aa='N'; codon[3][3][2].property='E';
	codon[3][3][3].aa='K'; codon[3][3][3].property='P'; codon[3][3][4].aa='K'; codon[3][3][4].property='P';
	codon[4][3][1].aa='D'; codon[4][3][1].property='N'; codon[4][3][2].aa='D'; codon[4][3][2].property='N';
	codon[4][3][3].aa='E'; codon[4][3][3].property='N'; codon[4][3][4].aa='E'; codon[4][3][4].property='N';
	codon[1][4][1].aa='C'; codon[1][4][1].property='E'; codon[1][4][2].aa='C'; codon[1][4][2].property='E';
	codon[1][4][3].aa='X'; codon[1][4][3].property='X'; codon[1][4][4].aa='W'; codon[1][4][4].property='E';
	codon[2][4][1].aa='R'; codon[2][4][1].property='P'; codon[2][4][2].aa='R'; codon[2][4][2].property='P';
	codon[2][4][3].aa='R'; codon[2][4][3].property='P'; codon[2][4][4].aa='R'; codon[2][4][4].property='P';
	codon[3][4][1].aa='S'; codon[3][4][1].property='E'; codon[3][4][2].aa='S'; codon[3][4][2].property='E';
	codon[3][4][3].aa='R'; codon[3][4][3].property='P'; codon[3][4][4].aa='R'; codon[3][4][4].property='P';
	codon[4][4][1].aa='G'; codon[4][4][1].property='E'; codon[4][4][2].aa='G'; codon[4][4][2].property='E';
	codon[4][4][3].aa='G'; codon[4][4][3].property='E'; codon[4][4][4].aa='G'; codon[4][4][4].property='E';



	if((fin=fopen(argv[1], "r"))==NULL) /****** sequence file ******/
	{
		printf("First sequence file not found\n");
		exit(1);
	}
	entnum=0;
	while(fgets(buffer, MAXBUF, fin)!=NULL)
	{
		if(!isspace(buffer[0]) && buffer[0]!='/')
		{
			for(aaa=15; isspace(buffer[aaa]); aaa++) { ; } /****** Skip the space ******/

			aab=0;
			for(aaa; buffer[aaa]!='\n'; aaa++) /****** Input the sequence ******/
			{
				if(!isspace(buffer[aaa])) /****** Omitting spaces, but gaps can be included ******/
				{
					if(buffer[aaa]=='U' || buffer[aaa]=='u' || buffer[aaa]=='T' || buffer[aaa]=='t')
					{ ++aab; }
					else if(buffer[aaa]=='C' || buffer[aaa]=='c')
					{ ++aab; }
					else if(buffer[aaa]=='A' || buffer[aaa]=='a')
					{ ++aab; }
					else if(buffer[aaa]=='G' || buffer[aaa]=='g')
					{ ++aab; }
					else if(buffer[aaa]=='-')
					{ printf("unknown nucleotide\n"); fflush(stdout); exit(1); }
				}
			}
			sitenum=aab; /* sitenum is the number of sites */
			++entnum; /* entnum is the number of entries */
		}
	}
	fclose(fin);



printf("%d sequences\n", entnum); fflush(stdout);
printf("%d nucleotide sites\n", sitenum); fflush(stdout);



        if((fin=fopen(argv[2], "r"))==NULL) /****** oyako file ******/
        {
                printf("Second oyako file not found\n"); fflush(stdout);
                exit(1);
        }
        for(aaa=0; aaa<MAXPAIR; aaa++)
        {
                pairs[aaa].name1[0]='\0';
                pairs[aaa].name2[0]='\0';
		pairs[aaa].blen=0.0;
        }
        npairs=0;
        while(fgets(buffer, MAXBUF, fin)!=NULL)
        {
                aab=0;
                for(aaa=0; !isspace(buffer[aaa]); aaa++)
                {
                        pairs[npairs].name1[aab]=buffer[aaa];
                        ++aab;
                }
                pairs[npairs].name1[aab]='\0';

                for(aaa; isspace(buffer[aaa]); aaa++) { ; }

                aab=0;
                for(aaa; !isspace(buffer[aaa]); aaa++)
                {
                        pairs[npairs].name2[aab]=buffer[aaa];
                        ++aab;
                }
                pairs[npairs].name2[aab]='\0';

                for(aaa; isspace(buffer[aaa]); aaa++) { ; }

		aab=0;
		for(aaa; !isspace(buffer[aaa]); aaa++)
		{
			blena[aab]=buffer[aaa];
			++aab;
		}
		blena[aab]='\0';

		pairs[npairs].blen=atof(blena);

                ++npairs;
        }
        fclose(fin);



printf("%d branches\n", npairs); fflush(stdout);



	for(aaa=1; aaa<=2000; aaa++)
	{
		accept[aaa][0]='\0';
	}
	if((fin=fopen(argv[3], "r"))==NULL) /****** observed amino acid file ******/
	{
		printf("The third observed amino acids file not found.\n"); fflush(stdout);
		exit(1);
	}
	aab=0;
	while(fgets(buffer, MAXBUF, fin)!=NULL)
	{
		++aab; /* codon site number starts from 1 */
		aac=0; /* acceptable amino acid number */
		for(aaa=0; !isspace(buffer[aaa]); aaa++) { ; }
		for(aaa; buffer[aaa]!='\n'; aaa++)
		{
			if(isalpha(buffer[aaa]))
			{
				accept[aab][aac]=buffer[aaa];
				++aac;
			}
			accept[aab][aac]='\0';
		}
	}
	fclose(fin);






for(aaa=1; aaa<=(sitenum/3); aaa++)
{
printf("%d", aaa); fflush(stdout);
for(aab=0; isalpha(accept[aaa][aab]); aab++)
{
printf(" %c", accept[aaa][aab]); fflush(stdout);
}
printf("\n"); fflush(stdout);
}








	for(aaa=0; aaa<sitenum; aaa=aaa+3) /***** for each codon site ******/
	{



		chgnumc1=0.0;
		chgnumc2=0.0;
		chgnumc3=0.0;
		chgnumaa=0.0;



        	for(saa=0; saa<MAXPAIR; saa++)
        	{
			for(baa=0; baa<5; baa++)
			{
               			pairs[saa].c1[baa]='\0';
               			pairs[saa].c2[baa]='\0';
               			pairs[saa].c3[baa]='\0';
			}
			for(baa=0; baa<21; baa++)
			{
				pairs[saa].aa[baa]='\0';
			}
		}

		if((fin=fopen(argv[1], "r"))==NULL)
       		{
       		       	printf("The first sequence file was not found\n");
                	exit(1);
        	}
		while(fgets(buffer, MAXBUF, fin)!=NULL)
		{
			bab=0;
			for(baa=0; !isspace(buffer[baa]); baa++)
			{
				bufname[bab]=buffer[baa];
				++bab;
			}
			bufname[bab]='\0';

			seqatta='n';
			for(saa=0; saa<npairs; saa++)
			{
				if(strcmp(pairs[saa].name2, bufname)==0)
				{
					seqatta='y';

					for(baa; isspace(buffer[baa]); baa++) { ; }

					bab=0;
					for(baa; isalpha(buffer[baa]); baa++)
					{
						bufseq[bab]=buffer[baa];
						++bab;
					}
					bufseq[bab]='\0';

					if(     bufseq[aaa]=='T' || bufseq[aaa]=='t' || bufseq[aaa]=='U' || bufseq[aaa]=='u') { cod1[0]=1; }
					else if(bufseq[aaa]=='C' || bufseq[aaa]=='c')                                         { cod1[0]=2; }
					else if(bufseq[aaa]=='A' || bufseq[aaa]=='a')                                         { cod1[0]=3; }
					else if(bufseq[aaa]=='G' || bufseq[aaa]=='g')                                         { cod1[0]=4; }
					if(     bufseq[aaa+1]=='T' || bufseq[aaa+1]=='t' || bufseq[aaa+1]=='U' || bufseq[aaa+1]=='u') { cod1[1]=1; }
					else if(bufseq[aaa+1]=='C' || bufseq[aaa+1]=='c')                                             { cod1[1]=2; }
					else if(bufseq[aaa+1]=='A' || bufseq[aaa+1]=='a')                                             { cod1[1]=3; }
					else if(bufseq[aaa+1]=='G' || bufseq[aaa+1]=='g')                                             { cod1[1]=4; }
					if(     bufseq[aaa+2]=='T' || bufseq[aaa+2]=='t' || bufseq[aaa+2]=='U' || bufseq[aaa+2]=='u') { cod1[2]=1; }
					else if(bufseq[aaa+2]=='C' || bufseq[aaa+2]=='c')                                             { cod1[2]=2; }
					else if(bufseq[aaa+2]=='A' || bufseq[aaa+2]=='a')                                             { cod1[2]=3; }
					else if(bufseq[aaa+2]=='G' || bufseq[aaa+2]=='g')                                             { cod1[2]=4; }
					cod1[3]='\0';

					bufamino=codon[cod1[0]][cod1[1]][cod1[2]].aa;





/****** observed amino acids ******/
					aminoatta='n';
					for(maa=0; accept[(aaa/3)+1][maa]!='\0'; maa++)
					{
						if(bufamino==accept[(aaa/3)+1][maa])
						{
							aminoatta='y';
						}
					}

					if(aminoatta=='y')
					{

						if(     bufseq[aaa]=='T' || bufseq[aaa]=='t' || bufseq[aaa]=='U' || bufseq[aaa]=='u') { pairs[saa].c1[0]='t'; }
						else if(bufseq[aaa]=='C' || bufseq[aaa]=='c')                                         { pairs[saa].c1[0]='c'; }
						else if(bufseq[aaa]=='A' || bufseq[aaa]=='a')                                         { pairs[saa].c1[0]='a'; }
						else if(bufseq[aaa]=='G' || bufseq[aaa]=='g')                                         { pairs[saa].c1[0]='g'; }
						if(     bufseq[aaa+1]=='T' || bufseq[aaa+1]=='t' || bufseq[aaa+1]=='U' || bufseq[aaa+1]=='u') { pairs[saa].c2[0]='t'; }
						else if(bufseq[aaa+1]=='C' || bufseq[aaa+1]=='c')                                             { pairs[saa].c2[0]='c'; }
						else if(bufseq[aaa+1]=='A' || bufseq[aaa+1]=='a')                                             { pairs[saa].c2[0]='a'; }
						else if(bufseq[aaa+1]=='G' || bufseq[aaa+1]=='g')                                             { pairs[saa].c2[0]='g'; }
						if(     bufseq[aaa+2]=='T' || bufseq[aaa+2]=='t' || bufseq[aaa+2]=='U' || bufseq[aaa+2]=='u') { pairs[saa].c3[0]='t'; }
						else if(bufseq[aaa+2]=='C' || bufseq[aaa+2]=='c')                                             { pairs[saa].c3[0]='c'; }
						else if(bufseq[aaa+2]=='A' || bufseq[aaa+2]=='a')                                             { pairs[saa].c3[0]='a'; }
						else if(bufseq[aaa+2]=='G' || bufseq[aaa+2]=='g')                                             { pairs[saa].c3[0]='g'; }

						pairs[saa].aa[0]=codon[cod1[0]][cod1[1]][cod1[2]].aa;

					}







/******
printf("%s %c%c%c %c\n", pairs[saa].name2, pairs[saa].c1[0], pairs[saa].c2[0], pairs[saa].c3[0], pairs[saa].aa[0]); fflush(stdout);
******/

				}
			}
			if(seqatta=='n')
			{
				printf("sequence not found\n"); fflush(stdout);
				exit(1);
			}
		}
		fclose(fin);







		for(paa=0; paa<npairs; paa++)
		{







			for(kac=0; kac<5; kac++)
			{
				kodomo[0].c1[kac]='\0';
				kodomo[0].c2[kac]='\0';
				kodomo[0].c3[kac]='\0';
				kodomo[1].c1[kac]='\0';
				kodomo[1].c2[kac]='\0';
				kodomo[1].c3[kac]='\0';
			}
			for(kac=0; kac<21; kac++)
			{
				kodomo[0].aa[kac]='\0';
				kodomo[1].aa[kac]='\0';
			}

			konum=0;
			for(pab=0; pab<npairs; pab++)
			{
				if(strcmp(pairs[paa].name2, pairs[pab].name1)==0)
				{

/******
printf("ok\n"); fflush(stdout);
******/

					if(konum==2)
					{
						printf("error1\n"); fflush(stdout);
						exit(1);
					}

					for(kac=0; kac<5; kac++)
					{
						kodomo[konum].c1[kac]=pairs[pab].c1[kac];
						kodomo[konum].c2[kac]=pairs[pab].c2[kac];
						kodomo[konum].c3[kac]=pairs[pab].c3[kac];
					}
					for(kac=0; kac<21; kac++)
					{
						kodomo[konum].aa[kac]=pairs[pab].aa[kac];
					}

					++konum;
				}
			}

			if(konum==2)
			{
/******
printf("%c%c%c %c%c%c\n", kodomo[0].c1[0], kodomo[0].c2[0], kodomo[0].c3[0], kodomo[1].c1[0], kodomo[1].c2[0], kodomo[1].c3[0]); fflush(stdout);
******/
/****** c1 ******/
				onajiattac1='n';
				for(kaa=0; kaa<5; kaa++) 
				{ for(kab=0; kab<5; kab++) 
					{ if(isalpha(kodomo[0].c1[kaa]) && isalpha(kodomo[1].c1[kab]) && (kodomo[0].c1[kaa]==kodomo[1].c1[kab]))
						{ onajiattac1='y'; } } }
				if(onajiattac1=='y')
				{
					numc1=0;
					for(kaa=0; kaa<5; kaa++)
					{ for(kab=0; kab<5; kab++)
						{ if(isalpha(kodomo[0].c1[kaa]) && isalpha(kodomo[1].c1[kab]) && (kodomo[0].c1[kaa]==kodomo[1].c1[kab]))
							{ pairs[paa].c1[numc1]=kodomo[0].c1[kaa]; ++numc1; } } }
				}
				else if(onajiattac1=='n')
				{
					numc1=0;
					for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c1[kaa])) { pairs[paa].c1[numc1]=kodomo[0].c1[kaa]; ++numc1; } }
					for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c1[kab])) { pairs[paa].c1[numc1]=kodomo[1].c1[kab]; ++numc1; } }
					if(isalpha(kodomo[0].c1[0]) && isalpha(kodomo[1].c1[0])) { chgnumc1=chgnumc1+1.0; }
				}
/****** c2 ******/
				onajiattac2='n';
				for(kaa=0; kaa<5; kaa++)
				{ for(kab=0; kab<5; kab++)
					{ if(isalpha(kodomo[0].c2[kaa]) && isalpha(kodomo[1].c2[kab]) && (kodomo[0].c2[kaa]==kodomo[1].c2[kab]))
						{ onajiattac2='y'; } } }
				if(onajiattac2=='y')
				{
					numc2=0;
					for(kaa=0; kaa<5; kaa++)
					{ for(kab=0; kab<5; kab++)
						{ if(isalpha(kodomo[0].c2[kaa]) && isalpha(kodomo[1].c2[kab]) && (kodomo[0].c2[kaa]==kodomo[1].c2[kab]))
							{ pairs[paa].c2[numc2]=kodomo[0].c2[kaa]; ++numc2; } } }
				}
				else if(onajiattac2=='n')
				{
					numc2=0;
					for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c2[kaa])) { pairs[paa].c2[numc2]=kodomo[0].c2[kaa]; ++numc2; } }
					for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c2[kab])) { pairs[paa].c2[numc2]=kodomo[1].c2[kab]; ++numc2; } }
					if(isalpha(kodomo[0].c2[0]) && isalpha(kodomo[1].c2[0])) { chgnumc2=chgnumc2+1.0; }
				}

/****** c3 ******/
				onajiattac3='n';
				for(kaa=0; kaa<5; kaa++)
				{ for(kab=0; kab<5; kab++)
					{ if(isalpha(kodomo[0].c3[kaa]) && isalpha(kodomo[1].c3[kab]) && (kodomo[0].c3[kaa]==kodomo[1].c3[kab]))
						{ onajiattac3='y'; } } }
				if(onajiattac3=='y')
				{
					numc3=0;
					for(kaa=0; kaa<5; kaa++)
					{ for(kab=0; kab<5; kab++)
						{ if(isalpha(kodomo[0].c3[kaa]) && isalpha(kodomo[1].c3[kab]) && (kodomo[0].c3[kaa]==kodomo[1].c3[kab]))
							{ pairs[paa].c3[numc3]=kodomo[0].c3[kaa]; ++numc3; } } }
				}
				else if(onajiattac3=='n')
				{
					numc3=0;
					for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c3[kaa])) { pairs[paa].c3[numc3]=kodomo[0].c3[kaa]; ++numc3; } }
					for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c3[kab])) { pairs[paa].c3[numc3]=kodomo[1].c3[kab]; ++numc3; } }
					if(isalpha(kodomo[0].c3[0]) && isalpha(kodomo[1].c3[0])) { chgnumc3=chgnumc3+1.0; }
				}
/****** aa ******/
				onajiattaaa='n';
				for(kaa=0; kaa<21; kaa++)
				{ for(kab=0; kab<21; kab++)
					{ if(isalpha(kodomo[0].aa[kaa]) && isalpha(kodomo[1].aa[kab]) && (kodomo[0].aa[kaa]==kodomo[1].aa[kab]))
						{ onajiattaaa='y'; } } }
				if(onajiattaaa=='y')
				{
					numaa=0;
					for(kaa=0; kaa<21; kaa++)
					{ for(kab=0; kab<21; kab++)
						{ if(isalpha(kodomo[0].aa[kaa]) && isalpha(kodomo[1].aa[kab]) && (kodomo[0].aa[kaa]==kodomo[1].aa[kab]))
							{ pairs[paa].aa[numaa]=kodomo[0].aa[kaa]; ++numaa; } } }
				}
				else if(onajiattaaa=='n')
				{
					numaa=0;
					for(kaa=0; kaa<21; kaa++) { if(isalpha(kodomo[0].aa[kaa])) { pairs[paa].aa[numaa]=kodomo[0].aa[kaa]; ++numaa; } }
					for(kab=0; kab<21; kab++) { if(isalpha(kodomo[1].aa[kab])) { pairs[paa].aa[numaa]=kodomo[1].aa[kab]; ++numaa; } }
					if(isalpha(kodomo[0].aa[0]) && isalpha(kodomo[1].aa[0])) { chgnumaa=chgnumaa+1.0; }
				}

			}
			else if(konum==1)
			{
				numc1=0;
				for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c1[kaa])) { pairs[paa].c1[numc1]=kodomo[0].c1[kaa]; ++numc1; } }
				for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c1[kab])) { pairs[paa].c1[numc1]=kodomo[1].c1[kab]; ++numc1; } }
				numc2=0;
				for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c2[kaa])) { pairs[paa].c2[numc2]=kodomo[0].c2[kaa]; ++numc2; } }
				for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c2[kab])) { pairs[paa].c2[numc2]=kodomo[1].c2[kab]; ++numc2; } }
				numc3=0;
				for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c3[kaa])) { pairs[paa].c3[numc3]=kodomo[0].c3[kaa]; ++numc3; } }
				for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c3[kab])) { pairs[paa].c3[numc3]=kodomo[1].c3[kab]; ++numc3; } }
				numaa=0;
				for(kaa=0; kaa<21; kaa++) { if(isalpha(kodomo[0].aa[kaa])) { pairs[paa].aa[numaa]=kodomo[0].aa[kaa]; ++numaa; } }
				for(kab=0; kab<21; kab++) { if(isalpha(kodomo[1].aa[kab])) { pairs[paa].aa[numaa]=kodomo[1].aa[kab]; ++numaa; } }
			}
			else if(konum==0)
			{
				;
/******
				printf("error1\n"); fflush(stdout);
******/
			}
			else if(konum>=3)
			{
				printf("error2\n"); fflush(stdout); exit(1);
			}




		}













		for(kac=0; kac<5; kac++)
		{
			kodomo[0].c1[kac]='\0';
			kodomo[0].c2[kac]='\0';
			kodomo[0].c3[kac]='\0';
			kodomo[1].c1[kac]='\0';
			kodomo[1].c2[kac]='\0';
			kodomo[1].c3[kac]='\0';
		}
		for(kac=0; kac<21; kac++)
		{
			kodomo[0].aa[kac]='\0';
			kodomo[1].aa[kac]='\0';
		}



		for(kac=0; kac<5; kac++)
		{
			kodomo[0].c1[kac]=pairs[npairs-2].c1[kac];
			kodomo[0].c2[kac]=pairs[npairs-2].c2[kac];
			kodomo[0].c3[kac]=pairs[npairs-2].c3[kac];
		}
		for(kac=0; kac<21; kac++)
		{
			kodomo[0].aa[kac]=pairs[npairs-2].aa[kac];
		}



		for(kac=0; kac<5; kac++)
		{
			kodomo[1].c1[kac]=pairs[npairs-1].c1[kac];
			kodomo[1].c2[kac]=pairs[npairs-1].c2[kac];
			kodomo[1].c3[kac]=pairs[npairs-1].c3[kac];
		}
		for(kac=0; kac<21; kac++)
		{
			kodomo[1].aa[kac]=pairs[npairs-1].aa[kac];
		}




/****** c1 ******/
		onajiattac1='n';
		for(kaa=0; kaa<5; kaa++)
		{
			for(kab=0; kab<5; kab++)
			{
				if(isalpha(kodomo[0].c1[kaa]) && isalpha(kodomo[1].c1[kab]) && (kodomo[0].c1[kaa]==kodomo[1].c1[kab]))
				{
					onajiattac1='y';
				}
			}
		}
		if(onajiattac1=='y')
		{
			numc1=0;
			for(kaa=0; kaa<5; kaa++)
			{
				for(kab=0; kab<5; kab++)
				{
					if(isalpha(kodomo[0].c1[kaa]) && isalpha(kodomo[1].c1[kab]) && (kodomo[0].c1[kaa]==kodomo[1].c1[kab]))
					{
						pairs[paa].c1[numc1]=kodomo[0].c1[kaa];
						++numc1;
					}
				}
			}
		}
		else if(onajiattac1=='n')
		{
			numc1=0;
			for(kaa=0; kaa<5; kaa++)
			{
				if(isalpha(kodomo[0].c1[kaa]))
				{
					pairs[paa].c1[numc1]=kodomo[0].c1[kaa];
					++numc1;
				}
			}
			for(kab=0; kab<5; kab++)
			{
				if(isalpha(kodomo[1].c1[kab]))
				{
					pairs[paa].c1[numc1]=kodomo[1].c1[kab];
					++numc1;
				}
			}
			if(isalpha(kodomo[0].c1[0]) && isalpha(kodomo[1].c1[0]))
			{
				chgnumc1=chgnumc1+1.0;
			}
		}
/****** c2 ******/
		onajiattac2='n';
		for(kaa=0; kaa<5; kaa++)
		{ for(kab=0; kab<5; kab++)
			{ if(isalpha(kodomo[0].c2[kaa]) && isalpha(kodomo[1].c2[kab]) && (kodomo[0].c2[kaa]==kodomo[1].c2[kab]))
				{ onajiattac2='y'; } } }
		if(onajiattac2=='y')
		{
			numc2=0;
			for(kaa=0; kaa<5; kaa++)
			{ for(kab=0; kab<5; kab++)
				{ if(isalpha(kodomo[0].c2[kaa]) && isalpha(kodomo[1].c2[kab]) && (kodomo[0].c2[kaa]==kodomo[1].c2[kab]))
					{ pairs[paa].c2[numc2]=kodomo[0].c2[kaa]; ++numc2; } } }
		}
		else if(onajiattac2=='n')
		{
			numc2=0;
			for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c2[kaa])) { pairs[paa].c2[numc2]=kodomo[0].c2[kaa]; ++numc2; } }
			for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c2[kab])) { pairs[paa].c2[numc2]=kodomo[1].c2[kab]; ++numc2; } }
			if(isalpha(kodomo[0].c2[0]) && isalpha(kodomo[1].c2[0]))
			{
				chgnumc2=chgnumc2+1.0;
			}
		}

/****** c3 ******/
		onajiattac3='n';
		for(kaa=0; kaa<5; kaa++)
		{ for(kab=0; kab<5; kab++)
			{ if(isalpha(kodomo[0].c3[kaa]) && isalpha(kodomo[1].c3[kab]) && (kodomo[0].c3[kaa]==kodomo[1].c3[kab]))
				{ onajiattac3='y'; } } }
		if(onajiattac3=='y')
		{
			numc3=0;
			for(kaa=0; kaa<5; kaa++)
			{ for(kab=0; kab<5; kab++)
				{ if(isalpha(kodomo[0].c3[kaa]) && isalpha(kodomo[1].c3[kab]) && (kodomo[0].c3[kaa]==kodomo[1].c3[kab]))
					{ pairs[paa].c3[numc3]=kodomo[0].c3[kaa]; ++numc3; } } }
		}
		else if(onajiattac3=='n')
		{
			numc3=0;
			for(kaa=0; kaa<5; kaa++) { if(isalpha(kodomo[0].c3[kaa])) { pairs[paa].c3[numc3]=kodomo[0].c3[kaa]; ++numc3; } }
			for(kab=0; kab<5; kab++) { if(isalpha(kodomo[1].c3[kab])) { pairs[paa].c3[numc3]=kodomo[1].c3[kab]; ++numc3; } }
			if(isalpha(kodomo[0].c3[0]) && isalpha(kodomo[1].c3[0]))
			{
				chgnumc3=chgnumc3+1.0;
			}
		}
/****** aa ******/
		onajiattaaa='n';
		for(kaa=0; kaa<21; kaa++)
		{ for(kab=0; kab<21; kab++)
			{ if(isalpha(kodomo[0].aa[kaa]) && isalpha(kodomo[1].aa[kab]) && (kodomo[0].aa[kaa]==kodomo[1].aa[kab]))
				{ onajiattaaa='y'; } } }
		if(onajiattaaa=='y')
		{
			numaa=0;
			for(kaa=0; kaa<21; kaa++)
			{ for(kab=0; kab<21; kab++)
				{ if(isalpha(kodomo[0].aa[kaa]) && isalpha(kodomo[1].aa[kab]) && (kodomo[0].aa[kaa]==kodomo[1].aa[kab]))
					{ pairs[paa].aa[numaa]=kodomo[0].aa[kaa]; ++numaa; } } }
		}
		else if(onajiattaaa=='n')
		{
			numaa=0;
			for(kaa=0; kaa<21; kaa++) { if(isalpha(kodomo[0].aa[kaa])) { pairs[paa].aa[numaa]=kodomo[0].aa[kaa]; ++numaa; } }
			for(kab=0; kab<21; kab++) { if(isalpha(kodomo[1].aa[kab])) { pairs[paa].aa[numaa]=kodomo[1].aa[kab]; ++numaa; } }
			if(isalpha(kodomo[0].aa[0]) && isalpha(kodomo[1].aa[0]))
			{
				chgnumaa=chgnumaa+1.0;
			}
		}






/******
for(saa=0; saa<npairs; saa++)
{
printf("%s", pairs[saa].name2); fflush(stdout);
printf(" "); fflush(stdout);
for(sab=0; pairs[saa].c1[sab]!='\0'; sab++) { printf("%c", pairs[saa].c1[sab]); fflush(stdout); }
printf(" "); fflush(stdout);
for(sab=0; pairs[saa].c2[sab]!='\0'; sab++) { printf("%c", pairs[saa].c2[sab]); fflush(stdout); }
printf(" "); fflush(stdout);
for(sab=0; pairs[saa].c3[sab]!='\0'; sab++) { printf("%c", pairs[saa].c3[sab]); fflush(stdout); }
printf(" "); fflush(stdout);
for(sab=0; pairs[saa].aa[sab]!='\0'; sab++) { printf("%c", pairs[saa].aa[sab]); fflush(stdout); }
printf("\n"); fflush(stdout);
}
******/

printf("%d %f %f %f %f\n", (aaa/3)+1, chgnumc1, chgnumc2, chgnumc3, chgnumaa); fflush(stdout);



/******
exit(1);
******/



















                for(zaa=1; zaa<=4; zaa++)
                {
                        for(zab=1; zab<=4; zab++)
                        {
                                for(zac=1; zac<=4; zac++)
                                {
                                        flag='n';
                                        for(zad=0; accept[(aaa/3)+1][zad]!='\0'; zad++)
                                        {
                                                if(codon[zaa][zab][zac].aa==accept[(aaa/3)+1][zad])
                                                {
                                                        flag='y';
                                                }
                                        }
                                        if(flag=='n')
                                        {
                                                codon[zaa][zab][zac].property=codon[zaa][zab][zac].aa;
                                        }
                                        else if(flag=='y')
                                        {
                                                codon[zaa][zab][zac].property='B';
                                        }
                                }
                        }
                }
                codon[1][3][3].property='X';
                codon[1][3][4].property='X';
                codon[1][4][3].property='X';



                wsynsite=0.0;
                wconsite=0.0;
                tsynchge=0.0;
                tconchge=0.0;
                for(aab=0; aab<npairs; aab++) /****** for each pair ******/
                {


                        seq1[0]='\0';
                        seq2[0]='\0';

                        atta1='n';
                        if((fin=fopen(argv[1], "r"))==NULL)
                        {
                                printf("The first sequence file was not found\n");
                                exit(1);
                        }
                        while(fgets(buffer, MAXBUF, fin)!=NULL)
                        {
                                bab=0;
                                for(baa=0; !isspace(buffer[baa]); baa++)
                                {
                                        bufname[bab]=buffer[baa];
                                        ++bab;
                                }
                                bufname[bab]='\0';

                                if(strcmp(pairs[aab].name2, bufname)==0)
                                {
                                        atta1='y';

                                        for(aac=15; isspace(buffer[aac]); aac++) { ; } /****** Skip the space ******/
                        
                                        aad=0;
                                        for(aac; buffer[aac]!='\n'; aac++) /****** Input the sequence ******/
                                        {
                                                if(!isspace(buffer[aac])) /****** Omitting spaces, gaps can be included ******/
                                                {
                                                        if(buffer[aac]=='U' || buffer[aac]=='u' || buffer[aac]=='T' || buffer[aac]=='t')
                                                        { seq1[aad]='t'; seq2[aad]='t'; ++aad; }
                                                        else if(buffer[aac]=='C' || buffer[aac]=='c')
                                                        { seq1[aad]='c'; seq2[aad]='c'; ++aad; }
                                                        else if(buffer[aac]=='A' || buffer[aac]=='a')
                                                        { seq1[aad]='a'; seq2[aad]='a'; ++aad; }
                                                        else if(buffer[aac]=='G' || buffer[aac]=='g')
                                                        { seq1[aad]='g'; seq2[aad]='g'; ++aad; }
                                                        else if(buffer[aac]=='-')
                                                        { printf("unknown nucleotide\n"); fflush(stdout); exit(1); }
                                                }
                                        }
                                        seq1[aad]='\0';
                                        seq2[aad]='\0';
                                }
                        }
                        fclose(fin);

/******
                        if(atta1=='n')
                        {
                                printf("error3\n"); fflush(stdout);
                                exit(1);
                        }
******/


/******
                        site=0.0;
                        diff=0.0;
                        for(aac=0; seq1[aac]!='\0'; aac++)
                        {
                                site=site+1.0;
                                if(seq1[aac]!=seq2[aac])
                                { diff=diff+1.0; }
                        }
                        branlen=diff/site;
******/






                        if(seq1[aaa]=='t')        { cod1[0]=1; }
                        else if(seq1[aaa]=='c')   { cod1[0]=2; }
                        else if(seq1[aaa]=='a')   { cod1[0]=3; }
                        else if(seq1[aaa]=='g')   { cod1[0]=4; }
                        if(seq1[aaa+1]=='t')      { cod1[1]=1; }
                        else if(seq1[aaa+1]=='c') { cod1[1]=2; }
                        else if(seq1[aaa+1]=='a') { cod1[1]=3; }
                        else if(seq1[aaa+1]=='g') { cod1[1]=4; }
                        if(seq1[aaa+2]=='t')      { cod1[2]=1; }
                        else if(seq1[aaa+2]=='c') { cod1[2]=2; }
                        else if(seq1[aaa+2]=='a') { cod1[2]=3; }
                        else if(seq1[aaa+2]=='g') { cod1[2]=4; }
                        cod1[3]='\0';

                        if(seq2[aaa]=='t')        { cod2[0]=1; }
                        else if(seq2[aaa]=='c')   { cod2[0]=2; }
                        else if(seq2[aaa]=='a')   { cod2[0]=3; }
                        else if(seq2[aaa]=='g')   { cod2[0]=4; }
                        if(seq2[aaa+1]=='t')      { cod2[1]=1; }
                        else if(seq2[aaa+1]=='c') { cod2[1]=2; }
                        else if(seq2[aaa+1]=='a') { cod2[1]=3; }
                        else if(seq2[aaa+1]=='g') { cod2[1]=4; }
                        if(seq2[aaa+2]=='t')      { cod2[2]=1; }
                        else if(seq2[aaa+2]=='c') { cod2[2]=2; }
                        else if(seq2[aaa+2]=='a') { cod2[2]=3; }
                        else if(seq2[aaa+2]=='g') { cod2[2]=4; }
                        cod2[3]='\0';



                        if(codon[cod1[0]][cod1[1]][cod1[2]].property=='B'
                        && codon[cod2[0]][cod2[1]][cod2[2]].property=='B')
                        {
                                funccod(cod1, cod2, accept[(aaa/3)+1], codres);

                                wsynsite=wsynsite+codres[0];
                                wconsite=wconsite+codres[4];
                                tsynchge=tsynchge+codres[2];
                                tconchge=tconchge+codres[6];
                        }











                }
                printf("%d     1.00000     %f     %f     1.00000     %f     %f     1.00000     1.00000\n",
                (aaa/3)+1, (wsynsite*3.0)/(wsynsite+wconsite), (chgnumc1+chgnumc2+chgnumc3-chgnumaa), (wconsite*3.0)/(wsynsite+wconsite), chgnumaa); fflush(stdout);







/******
printf("%d %f %f %f %f\n", (aaa/3)+1, chgnumc1, chgnumc2, chgnumc3, chgnumaa); fflush(stdout);
******/

















        }
}




































/******
        for(aaa=0; aaa<sitenum; aaa=aaa+3)
        {
        	for(aaa=0; aaa<MAXPAIR; aaa++)
        	{
               		pairs[aaa].c1[0]='\0';
               		pairs[aaa].c2[0]='\0';
               		pairs[aaa].c3[0]='\0';
			pairs[aaa].aa[0]='\0';
        	}
		brannum=0;
	        for(zaa=1; zaa<=4; zaa++)
       		{
                	for(zab=1; zab<=4; zab++)
                	{
                        	for(zac=1; zac<=4; zac++)
                        	{
                                	flag='n';
                                	for(zad=0; accept[(aaa/3)+1][zad]!='\0'; zad++)
                                	{
                                        	if(codon[zaa][zab][zac].aa==accept[(aaa/3)+1][zad])
                                        	{
                                               		flag='y';
                                        	}
                                	}
                                	if(flag=='n')
                                	{
                                        	codon[zaa][zab][zac].property=codon[zaa][zab][zac].aa;
                                	}
                                	else if(flag=='y')
                                	{
                                        	codon[zaa][zab][zac].property='B';
                                	}
                        	}
                	}
        	}
        	codon[1][3][3].property='X';
       		codon[1][3][4].property='X';
        	codon[1][4][3].property='X';
******/
/******
                wsynsite=0.0;
                wnonsite=0.0;
                tsynchge=0.0;
                tnonchge=0.0;
******/
/******
                wsynsite=0.0;
                wconsite=0.0;
                tsynchge=0.0;
                tconchge=0.0;
                for(aab=0; aab<npairs; aab++)
                {
			branlen=pairs[aab].blen;
******/
/******
			if(branlen<(1.0/1698.0))
			{
				goto owari1;
			}
			if(pairs[aab].name1[0]!='s' || pairs[aab].name2[0]!='s')
			{
				goto owari1;
			}
******/
/******
			++brannum;
			seq1[0]='\0';
			seq2[0]='\0';
			atta1='n';
			atta2='n';
		        if((fin=fopen(argv[1], "r"))==NULL)
       			{
       		        	printf("The first file was not found\n");
                		exit(1);
        		}
			while(fgets(buffer, MAXBUF, fin)!=NULL)
			{
				bab=0;
				for(baa=0; !isspace(buffer[baa]); baa++)
				{
					bufname[bab]=buffer[baa];
					++bab;
				}
				bufname[bab]='\0';
				if(strcmp(pairs[aab].name1, bufname)==0)
				{
					atta1='y';
                        		for(aac=15; isspace(buffer[aac]); aac++) { ; }
                        		aad=0;
                        		for(aac; buffer[aac]!='\n'; aac++)
                        		{
                                		if(!isspace(buffer[aac]))
                                		{
                                        		if(buffer[aac]=='U' || buffer[aac]=='u' || buffer[aac]=='T' || buffer[aac]=='t')
                                        		{ seq1[aad]='t'; ++aad; }
                                        		else if(buffer[aac]=='C' || buffer[aac]=='c')
                                        		{ seq1[aad]='c'; ++aad; }
                                        		else if(buffer[aac]=='A' || buffer[aac]=='a')
                                        		{ seq1[aad]='a'; ++aad; }
                                        		else if(buffer[aac]=='G' || buffer[aac]=='g')
                                        		{ seq1[aad]='g'; ++aad; }
                                        		else if(buffer[aac]=='-')
                                        		{ printf("unknown nucleotide\n"); fflush(stdout); exit(1); }
                                		}
                        		}
                        		seq1[aad]='\0';
				}
				if(strcmp(pairs[aab].name2, bufname)==0)
				{
					atta2='y';
                        		for(aac=15; isspace(buffer[aac]); aac++) { ; }
                        		aad=0;
                        		for(aac; buffer[aac]!='\n'; aac++)
                        		{
                                		if(!isspace(buffer[aac]))
                                		{
                                        		if(buffer[aac]=='U' || buffer[aac]=='u' || buffer[aac]=='T' || buffer[aac]=='t')
                                        		{ seq2[aad]='t'; ++aad; }
                                        		else if(buffer[aac]=='C' || buffer[aac]=='c')
                                        		{ seq2[aad]='c'; ++aad; }
                                        		else if(buffer[aac]=='A' || buffer[aac]=='a')
                                        		{ seq2[aad]='a'; ++aad; }
                                        		else if(buffer[aac]=='G' || buffer[aac]=='g')
                                        		{ seq2[aad]='g'; ++aad; }
                                        		else if(buffer[aac]=='-')
                                        		{ printf("unknown nucleotide\n"); fflush(stdout); exit(1); }
                                		}
                        		}
                        		seq2[aad]='\0';
				}
			}
			fclose(fin);
			if(atta1=='n' || atta2=='n')
			{
				printf("error\n"); fflush(stdout);
				exit(1);
			}
******/
/******
			site=0.0;
			diff=0.0;
			for(aac=0; seq1[aac]!='\0'; aac++)
			{
				site=site+1.0;
				if(seq1[aac]!=seq2[aac])
				{ diff=diff+1.0; }
			}
			branlen=diff/site;
******/
/******
                        if(seq1[aaa]=='t') { cod1[0]=1; }
                        else if(seq1[aaa]=='c') { cod1[0]=2; }
                        else if(seq1[aaa]=='a') { cod1[0]=3; }
                        else if(seq1[aaa]=='g') { cod1[0]=4; }
                        if(seq1[aaa+1]=='t') { cod1[1]=1; }
                        else if(seq1[aaa+1]=='c') { cod1[1]=2; }
                        else if(seq1[aaa+1]=='a') { cod1[1]=3; }
                        else if(seq1[aaa+1]=='g') { cod1[1]=4; }
                        if(seq1[aaa+2]=='t') { cod1[2]=1; }
                        else if(seq1[aaa+2]=='c') { cod1[2]=2; }
                        else if(seq1[aaa+2]=='a') { cod1[2]=3; }
                        else if(seq1[aaa+2]=='g') { cod1[2]=4; }
                        cod1[3]='\0';
                        if(seq2[aaa]=='t') { cod2[0]=1; }
                        else if(seq2[aaa]=='c') { cod2[0]=2; }
                        else if(seq2[aaa]=='a') { cod2[0]=3; }
                        else if(seq2[aaa]=='g') { cod2[0]=4; }
                        if(seq2[aaa+1]=='t') { cod2[1]=1; }
                        else if(seq2[aaa+1]=='c') { cod2[1]=2; }
                        else if(seq2[aaa+1]=='a') { cod2[1]=3; }
                        else if(seq2[aaa+1]=='g') { cod2[1]=4; }
                        if(seq2[aaa+2]=='t') { cod2[2]=1; }
                        else if(seq2[aaa+2]=='c') { cod2[2]=2; }
                        else if(seq2[aaa+2]=='a') { cod2[2]=3; }
                        else if(seq2[aaa+2]=='g') { cod2[2]=4; }
                        cod2[3]='\0';
******/
/******
			if(codon[cod1[0]][cod1[1]][cod1[2]].property=='B')
******/
/******
			if(codon[cod1[0]][cod1[1]][cod1[2]].property=='B'
			&& codon[cod2[0]][cod2[1]][cod2[2]].property=='B')
                        {
				funccod(cod1, cod2, accept[(aaa/3)+1], codres);

				wsynsite=wsynsite+codres[0]*branlen;
				wconsite=wconsite+codres[4]*branlen;
				tsynchge=tsynchge+codres[2];
				tconchge=tconchge+codres[6];
			}
******/
/******
                        wsynsite=wsynsite+codres[0]*branlen;
                        wnonsite=wnonsite+codres[1]*branlen;
                        tsynchge=tsynchge+codres[2];
                        tnonchge=tnonchge+codres[3];
******/
/******
			owari1: ;
                }
                printf("%d     1.00000     %f     %f     1.00000     %f     %f     1.00000     1.00000   %d\n",
******/
/******
                aaa/3+1, (wsynsite*3.0)/(wsynsite+wnonsite), tsynchge, (wnonsite*3.0)/(wsynsite+wnonsite), tnonchge); fflush(stdout);
******/
/******
                aaa/3+1, (wsynsite*3.0)/(wsynsite+wconsite), tsynchge, (wconsite*3.0)/(wsynsite+wconsite), tconchge, brannum); fflush(stdout);
        }
}
******/
