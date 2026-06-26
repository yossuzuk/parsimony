#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "sbsstr.h"

main(int argc, char *argv[])
{

	char seq1[MAXSEQ];
	char seq2[MAXSEQ];

	struct codons{char aa; char property;};
	struct codons codon[5][5][5];

        struct pair {char name1[20]; char name2[20]; double blen;};
        struct pair pairs[MAXPAIR];

	FILE *fin;
	char buffer[MAXBUF];
/******
	char sequon[MAXSEQ];
	char bufnodea1[MAXDIG], bufnodea2[MAXDIG];
******/
	char bufname[20];
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


	if(argc!=5)
	{
		printf("sbsstr.x : a program of adaptsite-p eliminating structural constraints\n"); fflush(stdout);
		printf("use      : sbsstr.x pamlseqfile oyakofile acceptableaminoacids kappa\n"); fflush(stdout);
		printf("example  : sbsstr.x pamlseqfile oyakofile acceptable 6.0\n"); fflush(stdout);
		exit(1);
	}

/******
	for(aaa=0; aaa<MAXSEQ; aaa++)
	{
		sequon[aaa]='N';
	}
******/

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



	if((fin=fopen(argv[1], "r"))==NULL)
	{
		printf("First file not found\n");
		exit(1);
	}

	entnum=1;
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
			++entnum; /* entnum is the number of entries + 1 */
		}
	}
	fclose(fin);



printf("%d sequences\n", entnum-1); fflush(stdout);
printf("%d nucleotide sites\n", sitenum); fflush(stdout);



        if((fin=fopen(argv[2], "r"))==NULL)
        {
                printf("Second file not found\n"); fflush(stdout);
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






	for(aaa=1; aaa<=2000; aaa++)
	{
		accept[aaa][0]='\0';
	}
	if((fin=fopen(argv[3], "r"))==NULL)
	{
		printf("The third file not found.\n"); fflush(stdout);
		exit(1);
	}
	aab=0;
	while(fgets(buffer, MAXBUF, fin)!=NULL)
	{
		++aab; /* codon number starts from 1 */
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




        for(aaa=0; aaa<sitenum; aaa=aaa+3) /****** for each site ******/
        {



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



/******
                wsynsite=0.0;
                wnonsite=0.0;
                tsynchge=0.0;
                tnonchge=0.0;
******/
                wsynsite=0.0;
                wconsite=0.0;
                tsynchge=0.0;
                tconchge=0.0;
                for(aab=0; aab<npairs; aab++) /****** for each pair ******/
                {
			branlen=pairs[aab].blen;


/******/
			if(branlen<(1.0/1698.0))
			{
				goto owari1;
			}
			if(pairs[aab].name1[0]!='s' || pairs[aab].name2[0]!='s')
			{
				goto owari1;
			}
/******/




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

                        		for(aac=15; isspace(buffer[aac]); aac++) { ; } /****** Skip the space ******/
			
                        		aad=0;
                        		for(aac; buffer[aac]!='\n'; aac++) /****** Input the sequence ******/
                        		{
                                		if(!isspace(buffer[aac])) /****** Omitting spaces, gaps can be included ******/
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

                        		for(aac=15; isspace(buffer[aac]); aac++) { ; } /****** Skip the space ******/
			
                        		aad=0;
                        		for(aac; buffer[aac]!='\n'; aac++) /****** Input the sequence ******/
                        		{
                                		if(!isspace(buffer[aac])) /****** Omitting spaces, gaps can be included ******/
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


/******
			if(codon[cod1[0]][cod1[1]][cod1[2]].property=='B')
******/
			if(codon[cod1[0]][cod1[1]][cod1[2]].property=='B'
			&& codon[cod2[0]][cod2[1]][cod2[2]].property=='B')
                        {
				funccod(cod1, cod2, accept[(aaa/3)+1], codres);

				wsynsite=wsynsite+codres[0]*branlen;
				wconsite=wconsite+codres[4]*branlen;
				tsynchge=tsynchge+codres[2];
				tconchge=tconchge+codres[6];
			}



/*
                        wsynsite=wsynsite+codres[0]*branlen;
                        wnonsite=wnonsite+codres[1]*branlen;
                        tsynchge=tsynchge+codres[2];
                        tnonchge=tnonchge+codres[3];
*/



/******/
			owari1: ;
/******/





                }
                printf("%d     1.00000     %f     %f     1.00000     %f     %f     1.00000     1.00000   %d\n",
/*
                aaa/3+1, (wsynsite*3.0)/(wsynsite+wnonsite), tsynchge, (wnonsite*3.0)/(wsynsite+wnonsite), tnonchge); fflush(stdout);
*/
                aaa/3+1, (wsynsite*3.0)/(wsynsite+wconsite), tsynchge, (wconsite*3.0)/(wsynsite+wconsite), tconchge, brannum); fflush(stdout);
        }
}
