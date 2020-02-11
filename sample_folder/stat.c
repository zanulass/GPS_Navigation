/*------------------------------------------------------------
 * STAT.c - Position Statistics.
 *------------------------------------------------------------*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*------------------------------------------------------------
 * �萔�E�\���̂̒�`
 *------------------------------------------------------------*/
/* �_���^ */
typedef	int	bool;
#define	TRUE	1
#define	FALSE	0

/* WGS-84�萔 */
#define	PI		3.1415926535898		/* �~�����iIS-GPS-200�j */
#define	C		2.99792458e8		/* ����[m/s] */
#define	MUe		3.986005e14			/* �n���d�͒萔[m^3/s^2] */
#define	dOMEGAe	7.2921151467e-5		/* �n�����]�p���x[rad/s] */
#define	Re		6378137.0			/* �n�����a[m] */
#define	Fe		(1.0/298.257223563)	/* �n���̝G���� */

/* �p�x�̕ϊ� */
#define	rad_to_deg(rad)		((rad)/PI*180.0)
#define	deg_to_rad(deg)		((deg)/180.0*PI)
#define	rad_to_sc(rad)		((rad)/PI)
#define	sc_to_rad(sc)		((sc)*PI)

/* ��舵����s��̑傫�� */
#define	MAX_N				16		/* �ϑ��q�����̏�� */
#define	MAX_M				4		/* ���m���̍ő吔 */
#define	MAX_PRN				32		/* �q���ԍ��̏�� */

/* ���� */
#define	SECONDS_DAY			(3600L*24L)
#define	SECONDS_WEEK		(3600L*24L*7L)

/* ������\���\���� */
typedef struct {
	int		week;		/* �T�ԍ� */
	double	sec;		/* �T���߂���̌o�ߎ���[s] */
}	wtime;

/* �������W��\���\���� */
typedef struct {
	double	x;			/* X���W[m] */
	double	y;			/* Y���W[m] */
	double	z;			/* Z���W[m] */
}	posxyz;
#define	SQ(x)			((x)*(x))
#define	DIST(a,b)		sqrt(SQ(a.x-b.x)+SQ(a.y-b.y)+SQ(a.z-b.z))

/* �o�ܓx��\���\���� */
typedef struct {
	double	lat;		/* �ܓx[rad] */
	double	lon;		/* �o�x[rad] */
	double	hgt;		/* ���x�i�ȉ~�̍��j[m] */
}	posblh;

/* ENU���W��\���\���� */
typedef struct {
	double	e;			/* East����[m] */
	double	n;			/* North����[m] */
	double	u;			/* Up����[m] */
}	posenu;

/*------------------------------------------------------------
 * �t�@�C������̓Ǎ���
 *------------------------------------------------------------*/
/* �s�o�b�t�@ */
#define	LINEBUF_LEN			256
static char	linebuf[LINEBUF_LEN],fieldbuf[LINEBUF_LEN];
static int	linepos			=0;

/* �t�@�C������s�o�b�t�@�ɓǂݍ��� */
static bool read_line(FILE *fp)
{
	/* �s�o�b�t�@���N���A */
	linepos		=0;
	linebuf[0]	='\0';

	/* �t�@�C������ǂݍ��� */
	return (fgets(linebuf,LINEBUF_LEN,fp)!=NULL);
}

/* �s�o�b�t�@���當��������o���i�������w��A�܂���CSV�`���j */
static char *get_field(int width)
{
	int		i;
	bool	quotef=TRUE;

	/* width==0�Ȃ�CSV�`�� */
	if (width<1) {
		width	=LINEBUF_LEN-1;
		quotef	=FALSE;
	}

	/* �w�肳�ꂽ���������邢�̓R���}�܂ł�ǂݎ�� */
	for(i=0;i<width && linebuf[linepos+i]!='\0';i++) {
		if (linebuf[linepos+i]=='\n') break;
		if (!quotef && linebuf[linepos+i]==',') {
			linepos++;
			break;
		}
		if (linebuf[linepos+i]=='\"') {
			quotef=!quotef;
			continue;
		}
		fieldbuf[i]=linebuf[linepos+i];
	}
	linepos		+=i;		/* ���̈ʒu */
	fieldbuf[i]	='\0';

	return fieldbuf;
}

/*------------------------------------------------------------
 * ���W�ϊ�
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * xyz_to_blh() - �������W����o�ܓx�ւ̕ϊ�
 *------------------------------------------------------------
 *  posblh xyz_to_blh(pos); �o�ܓx
 *    posxyz pos; �������W�l
 *------------------------------------------------------------*/
posblh xyz_to_blh(posxyz pos)
{
	double	a,b,e,f,n,h,p,t,sint,cost;
	posblh	blh={0.0,0.0,-Re};

	/* ���_�̏ꍇ */
	if (pos.x==0.0 && pos.y==0.0 && pos.z==0.0) return blh;

	/* �ȉ~�̂̃p�����[�^ */
	f	=Fe;				/* �G���� */
	a	=Re;				/* �����a */
	b	=a*(1.0-f);			/* �Z���a */
	e	=sqrt(f*(2.0-f));	/* ���S�� */

	/* ���W�ϊ��̂��߂̃p�����[�^ */
	h	=a*a-b*b;
	p	=sqrt(pos.x*pos.x+pos.y*pos.y);
	t	=atan2(pos.z*a,p*b);
	sint=sin(t);
	cost=cos(t);

	/* �o�ܓx�ւ̕ϊ� */
	blh.lat =atan2(pos.z+h/b*sint*sint*sint,p-h/a*cost*cost*cost);
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	blh.lon	=atan2(pos.y,pos.x);
	blh.hgt	=(p/cos(blh.lat))-n;

	return blh;
}

/*------------------------------------------------------------
 * blh_to_xyz() - �o�ܓx���璼�����W�ւ̕ϊ�
 *------------------------------------------------------------
 *  posxyz blh_to_xyz(blh); �������W�l
 *    posblh blh; �o�ܓx
 *------------------------------------------------------------*/
posxyz blh_to_xyz(posblh blh)
{
	double	a,b,e,f,n;
	posxyz	pos;

	/* �ȉ~�̂̃p�����[�^ */
	f	=Fe;				/* �G���� */
	a	=Re;				/* �����a */
	b	=a*(1.0-f);			/* �Z���a */
	e	=sqrt(f*(2.0-f));	/* ���S�� */

	/* �������W�n�ւ̕ϊ� */
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	pos.x	=(n+blh.hgt)*cos(blh.lat)*cos(blh.lon);
	pos.y	=(n+blh.hgt)*cos(blh.lat)*sin(blh.lon);
	pos.z	=(n*(1.0-e*e)+blh.hgt)*sin(blh.lat);

	return pos;
}

/*------------------------------------------------------------
 * xyz_to_enu() - �������W����ENU���W�ւ̕ϊ�
 *------------------------------------------------------------
 *  posenu xyz_to_enu(pos,base); ENU���W�l
 *    posxyz pos;  �������W�l
 *    posxyz base; ��ʒu
 *------------------------------------------------------------*/
posenu xyz_to_enu(posxyz pos,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posenu	enu;

	/* ��ʒu����̑��Έʒu */
	pos.x	-=base.x;
	pos.y	-=base.y;
	pos.z	-=base.z;

	/* ��ʒu�̌o�ܓx */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* ���Έʒu����]������ENU���W�ɕϊ����� */
	enu.e	=-pos.x*s1+pos.y*c1;
	enu.n	=-pos.x*c1*s2-pos.y*s1*s2+pos.z*c2;
	enu.u	=pos.x*c1*c2+pos.y*s1*c2+pos.z*s2;

	return enu;
}

/*------------------------------------------------------------
 * enu_to_xyz() - ENU���W���璼�����W�ւ̕ϊ�
 *------------------------------------------------------------
 *  posxyz enu_to_xyz(enu,base); �������W�l
 *    posenu enu;  ENU���W�l
 *    posxyz base; ��ʒu
 *------------------------------------------------------------*/
posxyz enu_to_xyz(posenu enu,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posxyz	pos;

	/* ��ʒu�̌o�ܓx */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* ENU���W����]�����đ��Έʒu�ɕϊ����� */
	pos.x	=-enu.e*s1-enu.n*c1*s2+enu.u*c1*c2;
	pos.y	=enu.e*c1-enu.n*s1*s2+enu.u*s1*c2;
	pos.z	=enu.n*c2+enu.u*s2;

	/* ��ʒu�ɉ����� */
	pos.x	+=base.x;
	pos.y	+=base.y;
	pos.z	+=base.z;

	return pos;
}

/*------------------------------------------------------------
 * elevation() - �p�����߂�
 *------------------------------------------------------------
 *  double elevation(sat,usr); �p[rad]
 *    posxyz sat; �q���̈ʒu
 *    posxyz usr; ���[�U�ʒu
 *------------------------------------------------------------*/
double elevation(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU���W�ɕϊ����ċp�����߂� */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.u,sqrt(enu.e*enu.e+enu.n*enu.n));
}

/*------------------------------------------------------------
 * azimuth() - ���ʊp�����߂�
 *------------------------------------------------------------
 *  double azimuth(sat,usr); ���ʊp[rad]
 *    posxyz sat; �q���̈ʒu
 *    posxyz usr; ���[�U�ʒu
 *------------------------------------------------------------*/
double azimuth(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU���W�ɕϊ����ĕ��ʊp�����߂� */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.e,enu.n);
}

/*------------------------------------------------------------
 * main() - ���C��
 *------------------------------------------------------------*/
void main(int argc,char **argv)
{
	int		i;
	long	n=0L;
	double	d,data[6],data0[6],stat1[6],stat2[6],max[6],min[6];

	/* �R�}���h���C���������`�F�b�N */
	if (argc!=1) {
		printf("stat - position statistics.\n");
		printf("\n");
		printf("usage: stat <in >out\n");
		printf("\n");
		exit(0);
	}

	/* ������ */
	for(i=0;i<6;i++) {
		stat1[i]=0.0;
		stat2[i]=0.0;
	}

	/* �f�[�^�t�@�C����ǂݍ��� */
	while(read_line(stdin)) {
		/* �R�����g�s�͓ǂݔ�΂� */
		if (linebuf[0]=='#' || linebuf[0]=='%') continue;

		d		=atof(get_field(0));			/* ���� */
		data[0]	=atof(get_field(0));			/* X���W[m] */
		data[1]	=atof(get_field(0));			/* y���W[m] */
		data[2]	=atof(get_field(0));			/* Z���W[m] */
		/* �������� */
		data[3]	=sqrt(SQ(data[0])+SQ(data[1]));
		/* �������� */
		data[4]	=fabs(data[2]);
		/* �􉽋��� */
		data[5]	=sqrt(SQ(data[3])+SQ(data[2]));
		/* �ُ�l�Ȃ�X�L�b�v */
		if (d==0.0 && data[5]==0.0) continue;

		/* �a�̌v�Z */
		for(i=0;i<6;i++) {
			if (n<1) {
				data0[i]=data[i];
				data[i]	=0.0;
			} else {
				data[i]	-=data0[i];
			}
			stat1[i]+=data[i];
			stat2[i]+=SQ(data[i]);
			if (n<1 || data[i]>max[i]) max[i]=data[i];
			if (n<1 || data[i]<min[i]) min[i]=data[i];
		}
		n++;
	}

	/* �f�[�^�� */
	printf("NUM,%d\n",n);
	if (n<1) exit(0);

	/* �f�[�^���Ŋ��� */
	for(i=0;i<6;i++) {
		stat1[i]/=n;
		stat2[i]/=n;
	}

	/* ���ϒl */
	printf("AVG");
	for(i=0;i<6;i++) {
		printf(",%.4f",stat1[i]+data0[i]);
	}
	printf("\n");

	/* �W���΍� */
	printf("STD");
	for(i=0;i<6;i++) {
		printf(",%.4f",sqrt(stat2[i]-SQ(stat1[i])));
	}
	printf("\n");

	/* RMS */
	printf("RMS");
	for(i=0;i<6;i++) {
		printf(",%.4f",
			sqrt(SQ(data0[i])+2.0*data0[i]*stat1[i]+stat2[i]));
	}
	printf("\n");

	/* �ő�l */
	printf("MAX");
	for(i=0;i<6;i++) {
		printf(",%.4f",max[i]+data0[i]);
	}
	printf("\n");

	/* �ŏ��l */
	printf("MIN");
	for(i=0;i<6;i++) {
		printf(",%.4f",min[i]+data0[i]);
	}
	printf("\n");

	/* ���ψʒu */
	printf("OFF");
	for(i=0;i<3;i++) {
		printf(",%.4f",stat1[i]+data0[i]);
	}
	d=SQ(stat1[0]+data0[0])+SQ(stat1[1]+data0[1]);
	printf(",%.4f",sqrt(d));
	printf(",%.4f",fabs(stat1[2]+data0[2]));
	printf(",%.4f",sqrt(d+SQ(stat1[2]+data0[2])));
	printf("\n");

	exit(0);
}

