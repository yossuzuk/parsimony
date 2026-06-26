#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "sbsstr.h"

void funccod(char *codon3, char *codon4, char *accept, double *count)
{
	struct codons{int aa; int property; double synsite; double nonsite; double consite; double radsite;};
	struct codons codon[6][6][6];
	double possible1, possible2, possible3, pathway;
	double synmut1, synmut2, synmut3;
	double nonmut1, nonmut2, nonmut3;
	double conmut1, conmut2, conmut3;
	double radmut1, radmut2, radmut3;
	int aaa, aab, aac, aad;

	char flag;

	codon[1][1][1].aa='F'; codon[1][1][1].property='E';
	codon[1][1][2].aa='F'; codon[1][1][2].property='E';
	codon[1][1][3].aa='L'; codon[1][1][3].property='E';
	codon[1][1][4].aa='L'; codon[1][1][4].property='E';
	codon[2][1][1].aa='L'; codon[2][1][1].property='E';
	codon[2][1][2].aa='L'; codon[2][1][2].property='E';
	codon[2][1][3].aa='L'; codon[2][1][3].property='E';
	codon[2][1][4].aa='L'; codon[2][1][4].property='E';
	codon[3][1][1].aa='I'; codon[3][1][1].property='E';
	codon[3][1][2].aa='I'; codon[3][1][2].property='E';
	codon[3][1][3].aa='I'; codon[3][1][3].property='E';
	codon[3][1][4].aa='M'; codon[3][1][4].property='E';
	codon[4][1][1].aa='V'; codon[4][1][1].property='E';
	codon[4][1][2].aa='V'; codon[4][1][2].property='E';
	codon[4][1][3].aa='V'; codon[4][1][3].property='E';
	codon[4][1][4].aa='V'; codon[4][1][4].property='E';

	codon[1][2][1].aa='S'; codon[1][2][1].property='E';
	codon[1][2][2].aa='S'; codon[1][2][2].property='E';
	codon[1][2][3].aa='S'; codon[1][2][3].property='E';
	codon[1][2][4].aa='S'; codon[1][2][4].property='E';
	codon[2][2][1].aa='P'; codon[2][2][1].property='E';
	codon[2][2][2].aa='P'; codon[2][2][2].property='E';
	codon[2][2][3].aa='P'; codon[2][2][3].property='E';
	codon[2][2][4].aa='P'; codon[2][2][4].property='E';
	codon[3][2][1].aa='T'; codon[3][2][1].property='E';
	codon[3][2][2].aa='T'; codon[3][2][2].property='E';
	codon[3][2][3].aa='T'; codon[3][2][3].property='E';
	codon[3][2][4].aa='T'; codon[3][2][4].property='E';
	codon[4][2][1].aa='A'; codon[4][2][1].property='E';
	codon[4][2][2].aa='A'; codon[4][2][2].property='E';
	codon[4][2][3].aa='A'; codon[4][2][3].property='E';
	codon[4][2][4].aa='A'; codon[4][2][4].property='E';

	codon[1][3][1].aa='Y'; codon[1][3][1].property='E';
	codon[1][3][2].aa='Y'; codon[1][3][2].property='E';
	codon[1][3][3].aa='X'; codon[1][3][3].property='X';
	codon[1][3][4].aa='X'; codon[1][3][4].property='X';
	codon[2][3][1].aa='H'; codon[2][3][1].property='P';
	codon[2][3][2].aa='H'; codon[2][3][2].property='P';
	codon[2][3][3].aa='Q'; codon[2][3][3].property='E';
	codon[2][3][4].aa='Q'; codon[2][3][4].property='E';
	codon[3][3][1].aa='N'; codon[3][3][1].property='E';
	codon[3][3][2].aa='N'; codon[3][3][2].property='E';
	codon[3][3][3].aa='K'; codon[3][3][3].property='P';
	codon[3][3][4].aa='K'; codon[3][3][4].property='P';
	codon[4][3][1].aa='D'; codon[4][3][1].property='N';
	codon[4][3][2].aa='D'; codon[4][3][2].property='N';
	codon[4][3][3].aa='E'; codon[4][3][3].property='N';
	codon[4][3][4].aa='E'; codon[4][3][4].property='N';

	codon[1][4][1].aa='C'; codon[1][4][1].property='E';
	codon[1][4][2].aa='C'; codon[1][4][2].property='E';
	codon[1][4][3].aa='X'; codon[1][4][3].property='X';
	codon[1][4][4].aa='W'; codon[1][4][4].property='E';
	codon[2][4][1].aa='R'; codon[2][4][1].property='P';
	codon[2][4][2].aa='R'; codon[2][4][2].property='P';
	codon[2][4][3].aa='R'; codon[2][4][3].property='P';
	codon[2][4][4].aa='R'; codon[2][4][4].property='P';
	codon[3][4][1].aa='S'; codon[3][4][1].property='E';
	codon[3][4][2].aa='S'; codon[3][4][2].property='E';
	codon[3][4][3].aa='R'; codon[3][4][3].property='P';
	codon[3][4][4].aa='R'; codon[3][4][4].property='P';
	codon[4][4][1].aa='G'; codon[4][4][1].property='E';
	codon[4][4][2].aa='G'; codon[4][4][2].property='E';
	codon[4][4][3].aa='G'; codon[4][4][3].property='E';
	codon[4][4][4].aa='G'; codon[4][4][4].property='E';

	for(aaa=1; aaa<=4; aaa++)
	{
		for(aab=1; aab<=4; aab++)
		{
			for(aac=1; aac<=4; aac++)
			{
				flag='n';
				for(aad=0; accept[aad]!='\0'; aad++)
				{
					if(codon[aaa][aab][aac].aa==accept[aad])
					{
						flag='y';
					}
				}
				if(flag=='n')
				{
					codon[aaa][aab][aac].property=codon[aaa][aab][aac].aa;
				}
				else if(flag=='y')
				{
					codon[aaa][aab][aac].property='B';
				}
			}
		}
	}
	codon[1][3][3].property='X';
	codon[1][3][4].property='X';
	codon[1][4][3].property='X';





	for(aaa=1; aaa<=4; aaa++)
	{
		for(aab=1; aab<=4; aab++)
		{
			for(aac=1; aac<=4; aac++)
			{
				if(codon[aaa][aab][aac].aa!='X')
				{
					possible1=0.0; synmut1=0.0; nonmut1=0.0; conmut1=0.0; radmut1=0.0;
					for(aad=1; aad<=4; aad++) /* Mutation at the first position */
					{
						if(aaa!=aad)
						{
							if(codon[aad][aab][aac].aa!='X')
							{
								possible1=possible1+mu[aaa][aad];
								if(codon[aaa][aab][aac].aa==codon[aad][aab][aac].aa)
								{ synmut1=synmut1+mu[aaa][aad]; }
								else
								{
									nonmut1=nonmut1+mu[aaa][aad];
									if(codon[aaa][aab][aac].property==codon[aad][aab][aac].property)
									{ conmut1=conmut1+mu[aaa][aad]; }
									else
									{ radmut1=radmut1+mu[aaa][aad]; }
								}
							}
						}
					}
					possible2=0.0; synmut2=0.0; nonmut2=0.0; conmut2=0.0; radmut2=0.0;
					for(aad=1; aad<=4; aad++) /* Mutation at the second position */
					{
						if(aab!=aad)
						{
							if(codon[aaa][aad][aac].aa!='X')
							{
								possible2=possible2+mu[aab][aad];
								if(codon[aaa][aab][aac].aa==codon[aaa][aad][aac].aa)
								{ synmut2=synmut2+mu[aab][aad]; }
								else
								{
									nonmut2=nonmut2+mu[aab][aad];
									if(codon[aaa][aab][aac].property==codon[aaa][aad][aac].property)
									{ conmut2=conmut2+mu[aab][aad]; }
									else
									{ radmut2=radmut2+mu[aab][aad]; }
								}
							}
						}
					}
					possible3=0.0; synmut3=0.0; nonmut3=0.0; conmut3=0.0; radmut3=0.0;
					for(aad=1; aad<=4; aad++) /* Mutation at the third position */
					{
						if(aac!=aad)
						{
							if(codon[aaa][aab][aad].aa!='X')
							{
								possible3=possible3+mu[aac][aad];
								if(codon[aaa][aab][aac].aa==codon[aaa][aab][aad].aa)
								{ synmut3=synmut3+mu[aac][aad]; }
								else
								{
									nonmut3=nonmut3+mu[aac][aad];
									if(codon[aaa][aab][aac].property==codon[aaa][aab][aad].property)
									{ conmut3=conmut3+mu[aac][aad]; }
									else
									{ radmut3=radmut3+mu[aac][aad]; }
								}
							}
						}
					}
/*
					codon[aaa][aab][aac].synsite=(synmut1/possible1)+(synmut2/possible2)+(synmut3/possible3);
					codon[aaa][aab][aac].nonsite=(nonmut1/possible1)+(nonmut2/possible2)+(nonmut3/possible3);
					codon[aaa][aab][aac].consite=(conmut1/possible1)+(conmut2/possible2)+(conmut3/possible3);
					codon[aaa][aab][aac].radsite=(radmut1/possible1)+(radmut2/possible2)+(radmut3/possible3);
*/
					codon[aaa][aab][aac].synsite=synmut1+synmut2+synmut3;
					codon[aaa][aab][aac].nonsite=nonmut1+nonmut2+nonmut3;
					codon[aaa][aab][aac].consite=conmut1+conmut2+conmut3;
					codon[aaa][aab][aac].radsite=radmut1+radmut2+radmut3;
				}
			}
		}
	}







	/* any termination codon */
	if((codon3[0]==1 && codon3[1]==3 && codon3[2]==3) || (codon3[0]==1 && codon3[1]==3 && codon3[2]==4) || (codon3[0]==1 && codon3[1]==4 && codon3[2]==3) ||        (codon4[0]==1 && codon4[1]==3 && codon4[2]==3) || (codon4[0]==1 && codon4[1]==3 && codon4[2]==4) || (codon4[0]==1 && codon4[1]==4 && codon4[2]==3))
	{
		printf("error1\n"); fflush(stdout);
		exit(1);
	}

	pathway=0.0;
	count[0]=count[1]=count[2]=count[3]=count[4]=count[5]=count[6]=count[7]=0.0;

	/* identical */
	if(codon3[0]==codon4[0] && codon3[1]==codon4[1] && codon3[2]==codon4[2]) /* identical */
	{
		count[0]=codon[codon3[0]][codon3[1]][codon3[2]].synsite; /* synonymous site */
		count[1]=codon[codon3[0]][codon3[1]][codon3[2]].nonsite; /* nonsynonymous site */
		count[2]=0;                                              /* synonymous change */
		count[3]=0;                                              /* nonsynonymous change */
		count[4]=codon[codon3[0]][codon3[1]][codon3[2]].consite; /* conservative site */
		count[5]=codon[codon3[0]][codon3[1]][codon3[2]].radsite; /* radical site */
		count[6]=0;                                              /* conservative change */
		count[7]=0;                                              /* radical change */
		return;
	}

	/* one-position difference */
	else if((codon3[0]!=codon4[0] && codon3[1]==codon4[1] && codon3[2]==codon4[2]) || (codon3[0]==codon4[0] && codon3[1]!=codon4[1] && codon3[2]==codon4[2]) || (codon3[0]==codon4[0] && codon3[1]==codon4[1] && codon3[2]!=codon4[2])) /* one-position difference  */
	{
		count[0]=(codon[codon3[0]][codon3[1]][codon3[2]].synsite+codon[codon4[0]][codon4[1]][codon4[2]].synsite)/2.0; /* synonymous site */
		count[1]=(codon[codon3[0]][codon3[1]][codon3[2]].nonsite+codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/2.0; /* nonsynonymous site */
		count[4]=(codon[codon3[0]][codon3[1]][codon3[2]].consite+codon[codon4[0]][codon4[1]][codon4[2]].consite)/2.0; /* synonymous site */
		count[5]=(codon[codon3[0]][codon3[1]][codon3[2]].radsite+codon[codon4[0]][codon4[1]][codon4[2]].radsite)/2.0; /* nonsynonymous site */
		if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa) /* synonymous change */
		{
			count[2]=1.0;
			count[3]=0.0;
		}
		else
		{
			count[2]=0.0;
			count[3]=1.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property) /* conservative change */
			{
				count[6]=1.0;
				count[7]=0.0;
			}
			else
			{
				count[6]=0.0;
				count[7]=1.0;
			}
		}
		return;
	}

	/* first and second positions difference */
	else if(codon3[0]!=codon4[0] && codon3[1]!=codon4[1] && codon3[2]==codon4[2]) /* first and second positions difference */
	{
		if(codon[codon4[0]][codon3[1]][codon3[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon3[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon3[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon4[1]][codon3[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon4[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		count[0]=count[0]/pathway;
		count[1]=count[1]/pathway;
		count[2]=count[2]/pathway;
		count[3]=count[3]/pathway;
		count[4]=count[4]/pathway;
		count[5]=count[5]/pathway;
		count[6]=count[6]/pathway;
		count[7]=count[7]/pathway;
		return;
	}

	else if(codon3[0]!=codon4[0] && codon3[1]==codon4[1] && codon3[2]!=codon4[2]) /* first and third positions difference */
	{
		if(codon[codon4[0]][codon3[1]][codon3[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon3[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon3[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon3[1]][codon4[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon3[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon3[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		count[0]=count[0]/pathway;
		count[1]=count[1]/pathway;
		count[2]=count[2]/pathway;
		count[3]=count[3]/pathway;
		count[4]=count[4]/pathway;
		count[5]=count[5]/pathway;
		count[6]=count[6]/pathway;
		count[7]=count[7]/pathway;
		return;
	}

	else if(codon3[0]==codon4[0] && codon3[1]!=codon4[1] && codon3[2]!=codon4[2]) /* second and third positions difference */
	{
		if(codon[codon3[0]][codon4[1]][codon3[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon4[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon3[1]][codon4[2]].aa!='X')
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/3.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/3.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon3[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/3.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/3.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon3[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		count[0]=count[0]/pathway;
		count[1]=count[1]/pathway;
		count[2]=count[2]/pathway;
		count[3]=count[3]/pathway;
		count[4]=count[4]/pathway;
		count[5]=count[5]/pathway;
		count[6]=count[6]/pathway;
		count[7]=count[7]/pathway;
		return;
	}

	else if(codon3[0]!=codon4[0] && codon3[1]!=codon4[1] && codon3[2]!=codon4[2]) /* three positions difference */
	{
		if(codon[codon3[0]][codon4[1]][codon3[2]].aa!='X' && codon[codon3[0]][codon4[1]][codon4[2]].aa!='X') /* 2-3-1 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon4[1]][codon3[2]].consite
				          +codon[codon3[0]][codon4[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon3[2]].aa==codon[codon3[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon3[2]].property==codon[codon3[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon4[1]][codon3[2]].aa!='X' && codon[codon4[0]][codon4[1]][codon3[2]].aa!='X') /* 2-1-3 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon4[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon4[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon4[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon4[0]][codon3[1]][codon3[2]].aa!='X' && codon[codon4[0]][codon4[1]][codon3[2]].aa!='X') /* 1-2-3 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon3[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon3[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon3[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon4[1]][codon3[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon4[1]][codon3[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon4[0]][codon3[1]][codon3[2]].aa!='X' && codon[codon4[0]][codon3[1]][codon4[2]].aa!='X') /* 1-3-2 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon3[1]][codon3[2]].consite
				          +codon[codon4[0]][codon3[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon3[1]][codon3[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon3[1]][codon3[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon3[2]].aa==codon[codon4[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon3[2]].property==codon[codon4[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon3[1]][codon4[2]].aa!='X' && codon[codon4[0]][codon3[1]][codon4[2]].aa!='X') /* 3-1-2 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon3[1]][codon4[2]].consite
				          +codon[codon4[0]][codon3[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon3[1]][codon4[2]].aa==codon[codon4[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon4[2]].property==codon[codon4[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon4[0]][codon3[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon4[0]][codon3[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		if(codon[codon3[0]][codon3[1]][codon4[2]].aa!='X' && codon[codon3[0]][codon4[1]][codon4[2]].aa!='X') /* 3-2-1 */
		{
			pathway=pathway+1.0;
			count[0]=count[0]+(codon[codon3[0]][codon3[1]][codon3[2]].synsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].synsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].synsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].synsite)/4.0;
			count[1]=count[1]+(codon[codon3[0]][codon3[1]][codon3[2]].nonsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].nonsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].nonsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].nonsite)/4.0;
			count[4]=count[4]+(codon[codon3[0]][codon3[1]][codon3[2]].consite
				          +codon[codon3[0]][codon3[1]][codon4[2]].consite
				          +codon[codon3[0]][codon4[1]][codon4[2]].consite
				          +codon[codon4[0]][codon4[1]][codon4[2]].consite)/4.0;
			count[5]=count[5]+(codon[codon3[0]][codon3[1]][codon3[2]].radsite
				          +codon[codon3[0]][codon3[1]][codon4[2]].radsite
				          +codon[codon3[0]][codon4[1]][codon4[2]].radsite
				          +codon[codon4[0]][codon4[1]][codon4[2]].radsite)/4.0;
			if(codon[codon3[0]][codon3[1]][codon3[2]].aa==codon[codon3[0]][codon3[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon3[2]].property==codon[codon3[0]][codon3[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon3[1]][codon4[2]].aa==codon[codon3[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon3[1]][codon4[2]].property==codon[codon3[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
			if(codon[codon3[0]][codon4[1]][codon4[2]].aa==codon[codon4[0]][codon4[1]][codon4[2]].aa)
			{ count[2]=count[2]+1.0; }
			else
			{
				count[3]=count[3]+1.0;
				if(codon[codon3[0]][codon4[1]][codon4[2]].property==codon[codon4[0]][codon4[1]][codon4[2]].property)
				{ count[6]=count[6]+1.0; }
				else
				{ count[7]=count[7]+1.0; }
			}
		}
		count[0]=count[0]/pathway;
		count[1]=count[1]/pathway;
		count[2]=count[2]/pathway;
		count[3]=count[3]/pathway;
		count[4]=count[4]/pathway;
		count[5]=count[5]/pathway;
		count[6]=count[6]/pathway;
		count[7]=count[7]/pathway;
		return;
	}
}
