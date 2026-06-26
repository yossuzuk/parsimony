#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXDIG 1000
#define MAXBUF 1000
#define MAXSEQ 100000

int main(int argc, char *argv[])
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

	int aasite[1000];
	int totalaa;
	char aa1, aa2, aa3;
	char sequonatta;
	int baa;
	struct sequons{char seq[11][5]; int num[11];};
	struct sequons sequon[1000];

	if(argc!=2)
	{
		printf("ngssitefreq.x : a program to translate fastafile and find ngs\n"); fflush(stdout);
		printf("use           : ngssitefreq.x fastafile\n"); fflush(stdout);
		printf("example       : ngssitefreq.x ha.nuc.fas\n"); fflush(stdout);
		exit(1);
	}

	for(aaa=0; aaa<1000; aaa++)
	{
		aasite[aaa]=0;

		for(baa=0; baa<=10; baa++)
		{
			sequon[aaa].seq[baa][0]='\0'; sequon[aaa].num[baa]=0;
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

			for(aaa=1; (buffer[aaa]!='\n' && buffer[aaa]!='\0'); aaa++)
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

/******
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
******/
/******
			printf(">%s\n", bufname); fflush(stdout);
******/
			output=0; /* frame 1 */
			for(aaa=0; aaa<=sitenum-9; aaa=aaa+3)
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
				aa1=codon[cod[0]][cod[1]][cod[2]].aa;
				for(aab=3; aab<=5; aab++)
				{
					cod[3]='\0';
					if(bufseq[aaa+aab]=='T' || bufseq[aaa+aab]=='t' || bufseq[aaa+aab]=='U' || bufseq[aaa+aab]=='u') { cod[aab-3]=1; }
					else if(bufseq[aaa+aab]=='C' || bufseq[aaa+aab]=='c') { cod[aab-3]=2; }
					else if(bufseq[aaa+aab]=='A' || bufseq[aaa+aab]=='a') { cod[aab-3]=3; }
					else if(bufseq[aaa+aab]=='G' || bufseq[aaa+aab]=='g') { cod[aab-3]=4; }
					else { aminoacid='X'; goto frame1owari; }
				}
				aa2=codon[cod[0]][cod[1]][cod[2]].aa;
				for(aab=6; aab<=8; aab++)
				{
					cod[3]='\0';
					if(bufseq[aaa+aab]=='T' || bufseq[aaa+aab]=='t' || bufseq[aaa+aab]=='U' || bufseq[aaa+aab]=='u') { cod[aab-6]=1; }
					else if(bufseq[aaa+aab]=='C' || bufseq[aaa+aab]=='c') { cod[aab-6]=2; }
					else if(bufseq[aaa+aab]=='A' || bufseq[aaa+aab]=='a') { cod[aab-6]=3; }
					else if(bufseq[aaa+aab]=='G' || bufseq[aaa+aab]=='g') { cod[aab-6]=4; }
					else { aminoacid='X'; goto frame1owari; }
				}
				aa3=codon[cod[0]][cod[1]][cod[2]].aa;

				frame1owari: ;

				++output;

/******
printf("%c%c%c\n", aa1, aa2, aa3); fflush(stdout);
******/

				if(aa1=='N' && aa2!='P' && (aa3=='S' || aa3=='T'))
				{
					++aasite[output];

					sequonatta='N';
					for(baa=0; sequon[output].seq[baa][0]!='\0'; baa++)
					{
						if(sequon[output].seq[baa][0]==aa1
						&& sequon[output].seq[baa][1]==aa2
						&& sequon[output].seq[baa][2]==aa3)
						{
							++sequon[output].num[baa];
							sequonatta='Y';
							goto sequonowari;
						}
					}

					if(sequonatta=='N')
					{
						for(baa=0; baa<=10; baa++)
						{
							if(sequon[output].seq[baa][0]=='\0')
							{
								sequon[output].seq[baa][0]=aa1;
								sequon[output].seq[baa][1]=aa2;
								sequon[output].seq[baa][2]=aa3;

								++sequon[output].num[baa];
								goto sequonowari;
							}
						}
					}
				}
				sequonowari: ;
						











			}
		}
	}
	fclose(fin);

	for(aaa=1; aaa<=output; aaa++)
	{
		printf("%d %d\n", aaa, aasite[aaa]); fflush(stdout);
	}

	for(aaa=1; aaa<=output; aaa++)
	{
		if(1<=aaa && aaa<=17)
		{
			printf("%d %d", aaa-18, aasite[aaa]); fflush(stdout);

			if(aasite[aaa]==0)
			{
				printf("\n"); fflush(stdout);
			}
			if(aasite[aaa]>0)
			{
				for(baa=0; sequon[aaa].seq[baa][0]!='\0'; baa++)
				{
					printf(" %c%c%c %d", sequon[aaa].seq[baa][0], sequon[aaa].seq[baa][1], sequon[aaa].seq[baa][2], sequon[aaa].num[baa]); fflush(stdout);
				}
				printf("\n"); fflush(stdout);
			}
		}
		if(18<=aaa)
		{
			printf("%d %d", aaa-17, aasite[aaa]); fflush(stdout);

			if(aasite[aaa]==0)
			{
				printf("\n"); fflush(stdout);
			}
			if(aasite[aaa]>0)
			{
				for(baa=0; sequon[aaa].seq[baa][0]!='\0'; baa++)
				{
					printf(" %c%c%c %d", sequon[aaa].seq[baa][0], sequon[aaa].seq[baa][1], sequon[aaa].seq[baa][2], sequon[aaa].num[baa]); fflush(stdout);
				}
				printf("\n"); fflush(stdout);
			}
		}
	}
}
