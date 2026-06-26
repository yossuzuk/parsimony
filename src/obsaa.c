#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXDIG 1000
#define MAXBUF 1000
#define MAXSEQ 100000

main(int argc, char *argv[])
{
        FILE *fin;
        char buffer[MAXBUF];
	char bufname[MAXDIG];
        char bufseq[MAXSEQ];
        char bufrevseq[MAXSEQ];
        fpos_t pos;
        int aaa, aab;
	int output, downcut;
	int length;

	struct codons{char aa; char property;};
	struct codons codon[5][5][5];

	int sitenum, entnum, position;
	char aminoacid;
	char cod[4];

	int aasite[1000][30];
	int totalaa;
	double thresh;

	if(argc!=3)
	{
		printf("transseq.x : a program to translate fasta file\n"); fflush(stdout);
		printf("use        : transseq.x fastafile  thresh\n"); fflush(stdout);
		printf("example    : transseq.x ha.nuc.fas 0.05\n"); fflush(stdout);
		exit(1);
	}

	thresh=atof(argv[2]);

	for(aaa=0; aaa<1000; aaa++)
	{
		for(aab=0; aab<30; aab++)
		{
			aasite[aaa][aab]=0;
		}
	}

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
/******/
	codon[1][3][3].aa='X'; codon[1][3][3].property='X'; codon[1][3][4].aa='X'; codon[1][3][4].property='X';
/******/
/******
	codon[1][3][3].aa='.'; codon[1][3][3].property='.'; codon[1][3][4].aa='.'; codon[1][3][4].property='.';
******/
	codon[2][3][1].aa='H'; codon[2][3][1].property='P'; codon[2][3][2].aa='H'; codon[2][3][2].property='P';
	codon[2][3][3].aa='Q'; codon[2][3][3].property='E'; codon[2][3][4].aa='Q'; codon[2][3][4].property='E';
	codon[3][3][1].aa='N'; codon[3][3][1].property='E'; codon[3][3][2].aa='N'; codon[3][3][2].property='E';
	codon[3][3][3].aa='K'; codon[3][3][3].property='P'; codon[3][3][4].aa='K'; codon[3][3][4].property='P';
	codon[4][3][1].aa='D'; codon[4][3][1].property='N'; codon[4][3][2].aa='D'; codon[4][3][2].property='N';
	codon[4][3][3].aa='E'; codon[4][3][3].property='N'; codon[4][3][4].aa='E'; codon[4][3][4].property='N';
	codon[1][4][1].aa='C'; codon[1][4][1].property='E'; codon[1][4][2].aa='C'; codon[1][4][2].property='E';
/******/
	codon[1][4][3].aa='X'; codon[1][4][3].property='X'; codon[1][4][4].aa='W'; codon[1][4][4].property='E';
/******/
/******
	codon[1][4][3].aa='.'; codon[1][4][3].property='.'; codon[1][4][4].aa='W'; codon[1][4][4].property='E';
******/
	codon[2][4][1].aa='R'; codon[2][4][1].property='P'; codon[2][4][2].aa='R'; codon[2][4][2].property='P';
	codon[2][4][3].aa='R'; codon[2][4][3].property='P'; codon[2][4][4].aa='R'; codon[2][4][4].property='P';
	codon[3][4][1].aa='S'; codon[3][4][1].property='E'; codon[3][4][2].aa='S'; codon[3][4][2].property='E';
	codon[3][4][3].aa='R'; codon[3][4][3].property='P'; codon[3][4][4].aa='R'; codon[3][4][4].property='P';
	codon[4][4][1].aa='G'; codon[4][4][1].property='E'; codon[4][4][2].aa='G'; codon[4][4][2].property='E';
	codon[4][4][3].aa='G'; codon[4][4][3].property='E'; codon[4][4][4].aa='G'; codon[4][4][4].property='E';

	if((fin=fopen(argv[1], "r"))==NULL)
	{
		printf("The fasta file was not found\n");
		exit(1);
	}
	while(fgets(buffer, MAXBUF, fin)!=NULL)
	{
		sitenum=0;
		if(buffer[0]=='>')
		{
			fgetpos(fin, &pos);

			aab=0;
/*
			for(aaa=1; buffer[aaa]!='|'; aaa++)
*/
			for(aaa=1; !isspace(buffer[aaa]); aaa++)
			{
				bufname[aab]=buffer[aaa];
				++aab;
			}
			bufname[aab]='\0';

			aab=0;
			while(fgets(buffer, MAXBUF, fin)!=NULL)
			{
				if(buffer[0]!='>')
				{
					for(aaa=0; buffer[aaa]!='\n'; aaa++)
					{
						if(!isspace(buffer[aaa]))
						{
							bufseq[aab]=tolower(buffer[aaa]);
							++sitenum;
							++aab;
						}
					}
					bufseq[aab]='\0';
				}
				else if(buffer[0]=='>')
				{
					goto owari1;
				}
			}
			owari1: ;
       		        fsetpos(fin, &pos);

			aab=0;
			for(aaa=sitenum-1; aaa>=0; aaa--)
			{
				if(bufseq[aaa]=='T' || bufseq[aaa]=='t' || bufseq[aaa]=='U' || bufseq[aaa]=='u')
				{ bufrevseq[aab]='A'; ++aab; }
				else if(bufseq[aaa]=='C' || bufseq[aaa]=='c')
				{ bufrevseq[aab]='G'; ++aab; }
				else if(bufseq[aaa]=='A' || bufseq[aaa]=='a')
				{ bufrevseq[aab]='T'; ++aab; }
				else if(bufseq[aaa]=='G' || bufseq[aaa]=='g')
				{ bufrevseq[aab]='C'; ++aab; }
				else
				{ bufrevseq[aab]='X'; ++aab; }
			}
/******
			printf(">%s\n", bufname); fflush(stdout);
******/
			output=0; /* frame 1 */
			for(aaa=0; aaa<=sitenum-3; aaa=aaa+3)
			{
				for(aab=0; aab<=2; aab++)
				{
					cod[3]='\0';
					if(bufseq[aaa+aab]=='T' || bufseq[aaa+aab]=='t' || bufseq[aaa+aab]=='U' || bufseq[aaa+aab]=='u') { cod[aab]=1; }
					else if(bufseq[aaa+aab]=='C' || bufseq[aaa+aab]=='c') { cod[aab]=2; }
					else if(bufseq[aaa+aab]=='A' || bufseq[aaa+aab]=='a') { cod[aab]=3; }
					else if(bufseq[aaa+aab]=='G' || bufseq[aaa+aab]=='g') { cod[aab]=4; }
					else { aminoacid='X'; goto frame1owari; }
				}
				aminoacid=codon[cod[0]][cod[1]][cod[2]].aa;
				frame1owari: ;

				++output;

				if(aminoacid=='A') { ++aasite[output][1]; }
				if(aminoacid=='C') { ++aasite[output][2]; }
				if(aminoacid=='D') { ++aasite[output][3]; }
				if(aminoacid=='E') { ++aasite[output][4]; }
				if(aminoacid=='F') { ++aasite[output][5]; }
				if(aminoacid=='G') { ++aasite[output][6]; }
				if(aminoacid=='H') { ++aasite[output][7]; }
				if(aminoacid=='I') { ++aasite[output][8]; }
				if(aminoacid=='K') { ++aasite[output][9]; }
				if(aminoacid=='L') { ++aasite[output][10]; }
				if(aminoacid=='M') { ++aasite[output][11]; }
				if(aminoacid=='N') { ++aasite[output][12]; }
				if(aminoacid=='P') { ++aasite[output][13]; }
				if(aminoacid=='Q') { ++aasite[output][14]; }
				if(aminoacid=='R') { ++aasite[output][15]; }
				if(aminoacid=='S') { ++aasite[output][16]; }
				if(aminoacid=='T') { ++aasite[output][17]; }
				if(aminoacid=='V') { ++aasite[output][18]; }
				if(aminoacid=='W') { ++aasite[output][19]; }
				if(aminoacid=='Y') { ++aasite[output][20]; }
				if(aminoacid=='X') { ++aasite[output][21]; printf("error\n"); fflush(stdout); exit(1); }

/******
				if(output>0 && output%60==0) { printf("\n"); fflush(stdout); }
				printf("%c", aminoacid); fflush(stdout);
******/

			}
/******
			printf("\n"); fflush(stdout);
******/


		}
	}
	fclose(fin);

	for(aaa=1; aaa<=output; aaa++)
	{
		printf("%d", aaa); fflush(stdout);

		totalaa=0;
		for(aab=1; aab<=21; aab++)
		{
			totalaa=totalaa+aasite[aaa][aab];
		}
/******
		printf("%d ", totalaa); fflush(stdout);
******/

		for(aab=1; aab<=21; aab++)
		{
			if((((double)aasite[aaa][aab])/(double)totalaa)>thresh)
			{
/******
				if(aab==1) { printf("A %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==2) { printf("C %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==3) { printf("D %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==4) { printf("E %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==5) { printf("F %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==6) { printf("G %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==7) { printf("H %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==8) { printf("I %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==9) { printf("K %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==10) { printf("L %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==11) { printf("M %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==12) { printf("N %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==13) { printf("P %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==14) { printf("Q %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==15) { printf("R %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==16) { printf("S %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==17) { printf("T %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==18) { printf("V %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==19) { printf("W %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==20) { printf("Y %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
				if(aab==21) { printf("X %d ", aasite[aaa][aab]); totalaa=totalaa+aasite[aaa][aab]; } fflush(stdout);
******/
				if(aab==1) { printf(" A"); } fflush(stdout);
				if(aab==2) { printf(" C"); } fflush(stdout);
				if(aab==3) { printf(" D"); } fflush(stdout);
				if(aab==4) { printf(" E"); } fflush(stdout);
				if(aab==5) { printf(" F"); } fflush(stdout);
				if(aab==6) { printf(" G"); } fflush(stdout);
				if(aab==7) { printf(" H"); } fflush(stdout);
				if(aab==8) { printf(" I"); } fflush(stdout);
				if(aab==9) { printf(" K"); } fflush(stdout);
				if(aab==10) { printf(" L"); } fflush(stdout);
				if(aab==11) { printf(" M"); } fflush(stdout);
				if(aab==12) { printf(" N"); } fflush(stdout);
				if(aab==13) { printf(" P"); } fflush(stdout);
				if(aab==14) { printf(" Q"); } fflush(stdout);
				if(aab==15) { printf(" R"); } fflush(stdout);
				if(aab==16) { printf(" S"); } fflush(stdout);
				if(aab==17) { printf(" T"); } fflush(stdout);
				if(aab==18) { printf(" V"); } fflush(stdout);
				if(aab==19) { printf(" W"); } fflush(stdout);
				if(aab==20) { printf(" Y"); } fflush(stdout);
/******
				if(aab==21) { printf(" X"); printf("error\n"); exit(1); } fflush(stdout);
******/
			}
		}
		printf("\n"); fflush(stdout);
	}
}
