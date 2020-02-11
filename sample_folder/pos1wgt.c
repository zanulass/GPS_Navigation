/*------------------------------------------------------------
 * POS1WGT.c - Standalone Positioning.
 *------------------------------------------------------------
 *  Copyright (C) 2007 T. Sakai, ENRI <sakai@enri.go.jp>
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
 * ���s���̃p�����[�^
 *------------------------------------------------------------*/
/* �o�͌`�� */
static bool	output_nmea		=FALSE;
static bool	output_rel		=FALSE;
static bool	output_diff		=FALSE;
static posblh	rel_base;

/* �p�}�X�N[rad] */
static double	mask_angle	=0.0;

/* �c���`�F�b�N�̃X���b�V�����h[m] */
#define	DPSR_THRESHOLD		20.0

/* �W�I�C�h��[m] */
#define	GEOIDAL_HEIGHT		38.0

/* �d�ݕt���̃p�����[�^ */
#define	MIN_WEIGHT			0.001		/* �d�݂̍ŏ��l */
#define	VAR_ZENITH			(0.8*0.8)	/* �[�������̕��U[m^2] */

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
 * �����̕ϊ�
 *------------------------------------------------------------*/
/* �J�����_�l�̊J�n�N */
#define	TIME_T_BASE_YEAR	1970

/* 1980�N1��6��00:00:00�̃J�����_�l */
#define	TIME_T_ORIGIN		315964800L

/* mktime()�֐���GMT�Łigmtime()�֐��ɑΉ��j */
static time_t mktime2(struct tm *tm)
{
	int		i;
	long	days=0L;
	static int	days_month[]={
		31,28,31,30,31,30,31,31,30,31,30,31
	};

	/* �o�ߓ����𓾂� */
	for(i=TIME_T_BASE_YEAR;i<tm->tm_year+1900;i++) {
		days+=(i%4==0)?366:365;
	}
	for(i=1;i<tm->tm_mon+1;i++) {
		days+=days_month[i-1];
		if (i==2 && tm->tm_year%4==0) days++;
	}
	days+=tm->tm_mday-1;

	/* �J�����_�l��Ԃ� */
	return ((days*24+tm->tm_hour)*60+tm->tm_min)*60+tm->tm_sec;
}

/*------------------------------------------------------------
 * wtime_to_date() - �T�ԍ��E�b��������ւ̕ϊ�
 *------------------------------------------------------------
 *  struct tm wtime_to_date(wt); �����ւ̕ϊ�����
 *    wtime wt; �T�ԍ��E�b
 *------------------------------------------------------------*/
struct tm wtime_to_date(wtime wt)
{
	time_t	t;

	/* �������̌o�ߎ��Ԃ������āC�J�����_�l�𓾂� */
	t	=(long)wt.week*SECONDS_WEEK+TIME_T_ORIGIN
			+(long)((wt.sec>0.0)?wt.sec+0.5:wt.sec-0.5);

	return *gmtime(&t);
}

/*------------------------------------------------------------
 * date_to_wtime() - ��������T�ԍ��E�b�ւ̕ϊ�
 *------------------------------------------------------------
 *  wtime date_to_wtime(tmbuf); �T�ԍ��E�b�ւ̕ϊ�����
 *    struct tm tmbuf; �������w��
 *------------------------------------------------------------*/
wtime date_to_wtime(struct tm tmbuf)
{
	time_t	t;
	wtime	wt;

	/* �w�肳�ꂽ�����̃J�����_�l */
	t	=mktime2(&tmbuf);

	/* 1�T�Ԃ̕b���Ŋ��������Ɨ]�� */
	wt.week	=(t-TIME_T_ORIGIN)/SECONDS_WEEK;
	wt.sec	=(t-TIME_T_ORIGIN)%SECONDS_WEEK;

	return wt;
}

/*------------------------------------------------------------
 * �q�@���b�Z�[�W�̏���
 *------------------------------------------------------------*/
/* RINEX�t�@�C���̏�� */
#define	RINEX_POS_COMMENT		60
#define	RINEX_NAV_LINES			8
#define	RINEX_NAV_FIELDS_LINE	4

/* �L������G�t�F�����X�̍ő吔 */
#define	MAX_EPHMS				20

/* �G�t�F�����X�̗L������[h] */
#define	EPHEMERIS_EXPIRE		2.0

/* �G�t�F�����X���i�[���邽�߂̍\���� */
typedef	struct {
	int		week;		/* �T�ԍ� */
	double	data[RINEX_NAV_LINES*RINEX_NAV_FIELDS_LINE];
}	ephm_info;

/* �p�����[�^�ԍ����` */
enum ephm_para {
	EPHM_TOC,	EPHM_AF0,	EPHM_AF1,	EPHM_AF2,	/* line 1 */
	EPHM_IODE,	EPHM_Crs,	EPHM_d_n,	EPHM_M0,	/* line 2 */
	EPHM_Cuc,	EPHM_e,		EPHM_Cus,	EPHM_sqrtA,	/* line 3 */
	EPHM_TOE,	EPHM_Cic,	EPHM_OMEGA0,EPHM_Cis,	/* line 4 */
	EPHM_i0,	EPHM_Crc,	EPHM_omega,	EPHM_dOmega,/* line 5 */
	EPHM_di,	EPHM_CAonL2,EPHM_WEEK,	EPHM_L2P,	/* line 6 */
	EPHM_acc,	EPHM_health,EPHM_TGD,	EPHM_IODC,	/* line 7 */
	EPHM_TOT,	EPHM_FIT							/* line 8 */
};

/* �G�t�F�����X��ێ�����z�� */
static ephm_info	ephm_buf[MAX_PRN][MAX_EPHMS];
static int			ephm_count[MAX_PRN];
static int			current_ephm[MAX_PRN];
static int			current_week=-1;

/* �[�b�̏�� */
static int			leap_sec	=0;

/* �d���w�␳��� */
#define	IONO_PARAMETERS			4
static double	iono_alpha[IONO_PARAMETERS];
static double	iono_beta[IONO_PARAMETERS];

/* �����񁨎����̕ϊ� */
static double atof2(char *str)
{
	char	*p;

	/* 'D'��'E'�ɕϊ����� */
	for(p=str;*p!='\0';p++) {
		if (*p=='D' || *p=='d') *p='E';
	}

	/* �����ɕϊ����ĕԂ� */
	return atof(str);
}

/* �R�����g���𒲂ׂ� */
static bool is_comment(char *str)
{
	return (strncmp(linebuf+RINEX_POS_COMMENT,str,strlen(str))==0);
}

/*------------------------------------------------------------
 * read_RINEX_NAV() - RINEX�q�@�t�@�C����ǂݍ���
 *------------------------------------------------------------
 *  void read_RINEX_NAV(fp);
 *    FILE *fp; �ǂݍ��ރt�@�C��
 *------------------------------------------------------------*/
void read_RINEX_NAV(FILE *fp)
{
	int		i,j,n,prn,line;
	bool	noerr=FALSE;
	double	d,t;
	wtime	wt;
	struct tm	tmbuf;
	ephm_info	info;

	/* ������ */
	if (current_week<0) {
		for(i=0;i<IONO_PARAMETERS;i++) {
			iono_alpha[i]	=0.0;
			iono_beta[i]	=0.0;
		}
		for(i=0;i<MAX_PRN;i++) {
			ephm_count[i]	=0;
		}
	}

	/* �w�b�_���� */
	while(read_line(fp)) {
		if (is_comment("ION ALPHA")) {
			iono_alpha[0]	=atof2(get_field(14));
			iono_alpha[1]	=atof2(get_field(12));
			iono_alpha[2]	=atof2(get_field(12));
			iono_alpha[3]	=atof2(get_field(12));
		} else if (is_comment("ION BETA")) {
			iono_beta[0]	=atof2(get_field(14));
			iono_beta[1]	=atof2(get_field(12));
			iono_beta[2]	=atof2(get_field(12));
			iono_beta[3]	=atof2(get_field(12));
		} else if (is_comment("LEAP SECONDS")) {
			leap_sec		=atoi(get_field(6));
		} else if (is_comment("END OF HEADER")) break;
	}

	/* �{�̂�ǂ� */
	fprintf(stderr,"Reading RINEX NAV... ");
	while(read_line(fp)) {
		n=0;
		for(line=0;line<RINEX_NAV_LINES;line++) {
			/* �ŏ��̃f�[�^��ǂݍ��� */
			if (line==0) {
				/* �ŏ��̍s */
				prn				=atoi(get_field(2));
				tmbuf.tm_year	=atoi(get_field(3));
				if (tmbuf.tm_year<80) tmbuf.tm_year+=100;
				tmbuf.tm_mon	=atoi(get_field(3))-1;
				tmbuf.tm_mday	=atoi(get_field(3));
				tmbuf.tm_hour	=atoi(get_field(3));
				tmbuf.tm_min	=atoi(get_field(3));
				tmbuf.tm_sec	=0;
				wt				=date_to_wtime(tmbuf);
				wt.sec			+=atof(get_field(5));
				info.week		=wt.week;
				d				=wt.sec;
			} else {
				/* 2�s�ڈȍ~ */
				if (!read_line(fp)) goto ERROR;
				d	=atof2(get_field(22));
			}
			info.data[n++]=d;

			/* �c��̃f�[�^ */
			for(i=1;i<RINEX_NAV_FIELDS_LINE;i++) {
				d	=atof2(get_field(19));
				info.data[n++]=d;
			}
		}
		if (prn<1) continue;

		/* �T�ԍ��𑵂��� */
		t	=(info.week-info.data[EPHM_WEEK])*SECONDS_WEEK;
		info.week			=info.data[EPHM_WEEK];
		info.data[EPHM_TOC]	+=t;
		current_week		=info.week;

		/* ���łɓ������̂��Ȃ����ǂ��� */
		for(i=0;i<ephm_count[prn-1];i++) {
			/* �T�ԍ�����v������̂��Ώ� */
			if (ephm_buf[prn-1][i].week!=info.week) continue;

			/* IODC����v���Ă��邩 */
			if (ephm_buf[prn-1][i].data[EPHM_IODC]
					==info.data[EPHM_IODC]) {
				/* ���M�����̑������̂��c�� */
				if (info.data[EPHM_TOT]
						<ephm_buf[prn-1][i].data[EPHM_TOT]) {
					ephm_buf[prn-1][i]=info;
				}
				prn=0;
				break;
			}
		}
		if (prn<1) continue;

		/* �z��Ɋi�[����i���M�����̏����j */
		if (ephm_count[prn-1]>=MAX_EPHMS) {
			fprintf(stderr,"Too long NAV file.\n");
			goto NOERROR;
		}
		for(i=0;i<ephm_count[prn-1];i++) {
			t	=(ephm_buf[prn-1][i].week-info.week)*SECONDS_WEEK
					+ephm_buf[prn-1][i].data[EPHM_TOT];
			if (info.data[EPHM_TOT]<t) break;
		}
		for(j=ephm_count[prn-1];j>i;j--) {
			ephm_buf[prn-1][j]=ephm_buf[prn-1][j-1];
		}
		ephm_buf[prn-1][i]=info;
		ephm_count[prn-1]++;
	}
NOERROR:
	noerr=TRUE;

ERROR:
	/* �G�t�F�����X���Ȃ��ꍇ */
	if (current_week<0) {
		fprintf(stderr,"Error: No ephemeris information.\n");
		return;
	}

	/* �G�t�F�����X�̂���q���̐� */
	n=0; for(prn=1;prn<=MAX_PRN;prn++) {
		if (ephm_count[prn-1]>0) n++;
	}
	fprintf(stderr,"week %d: %d satellites\n",current_week,n);

	/* �r���Ńt�@�C�����I����Ă����ꍇ */
	if (!noerr) {
		fprintf(stderr,"Error: Unexpected EOF.\n");
	}
}

/*------------------------------------------------------------
 * set_ephemeris() - �G�t�F�����X���Z�b�g
 *------------------------------------------------------------
 *  bool set_ephemeris(prn,wt,iode); TRUE:�Z�b�g����
 *    int prn;  �q��PRN�ԍ�(1�`)
 *    wtime wt; �������w��
 *    int iode; IODE���w��/-1:�w��Ȃ�
 *------------------------------------------------------------*/
bool set_ephemeris(int prn,wtime wt,int iode)
{
	int		i;
	double	t,t0;

	/* �V�����G�t�F�����X����T�� */
	for(i=ephm_count[prn-1]-1;i>=0;i--) {
		/* �L���������ł��邱�� */
		t0	=(ephm_buf[prn-1][i].week-wt.week)*SECONDS_WEEK;
		t	=ephm_buf[prn-1][i].data[EPHM_TOC]+t0;
		if (wt.sec<t-EPHEMERIS_EXPIRE*3600.0-0.1 ||
				t+EPHEMERIS_EXPIRE*3600.0+0.1<wt.sec) continue;

		/* IODE���`�F�b�N */
		if (iode>=0) {
			if (ephm_buf[prn-1][i].data[EPHM_IODE]==iode) break;
			continue;
		}

		/* ��M���� */
		t	=ephm_buf[prn-1][i].data[EPHM_TOT]+t0;
		if (t<wt.sec+0.1) break;	/* wt���O�̎��� */
	}

	/* �J�����g���Ƃ��ăZ�b�g���� */
	current_ephm[prn-1]=i;
	return (i>=0);
}

/*------------------------------------------------------------
 * get_ephemeris() - �G�t�F�����X�̃p�����[�^�𓾂�
 *------------------------------------------------------------
 *  double get_ephemeris(prn,para); �p�����[�^�l
 *    int prn;  �q��PRN�ԍ�(1�`)
 *    int para; �p�����[�^�ԍ�(0�`)
 *------------------------------------------------------------
 *   ���O��set_ephemeris()�ɂ��G�t�F�����X���Z�b�g����Ă���
 * ����.
 *------------------------------------------------------------*/
double get_ephemeris(int prn,int para)
{
	if (ephm_count[prn-1]<1 || current_ephm[prn-1]<0) {
		fprintf(stderr,"Missing ephemeris: PRN=%d.\n",prn);
		exit(2);
	}

	return ephm_buf[prn-1][current_ephm[prn-1]].data[para];
}

/*------------------------------------------------------------
 * satellite_clock() - �q���N���b�N�덷���v�Z
 *------------------------------------------------------------
 *  double satellite_clock(prn,wt,psr); �q���N���b�N�덷[s]
 *    int prn;    �q��PRN�ԍ�(1�`)
 *    wtime wt;   �������w��
 *    double psr; �������̓`������[m]
 *------------------------------------------------------------
 *   ���O��set_ephemeris()�ɂ��G�t�F�����X���Z�b�g����Ă���
 * ����. ������get_ephemeris()���g�p����.
 *------------------------------------------------------------*/
double satellite_clock(int prn,wtime wt,double psr)
{
	int		i;
	double	tk,tk0,sqrtA,e,n,Ek,Mk,dt,tr;

	/* ���ԍ������߂� */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOE);
	tk	=tk0-psr/C;							/* �`�����Ԃ̕� */

	/* ���S�ߓ_�pEk[rad]�����߂� */
	sqrtA	=get_ephemeris(prn,EPHM_sqrtA);
	e	=get_ephemeris(prn,EPHM_e); 		/* ���S�� */
	n	=sqrt(MUe)/sqrtA/sqrtA/sqrtA+get_ephemeris(prn,EPHM_d_n);
	Mk	=get_ephemeris(prn,EPHM_M0)+n*tk;	/* ���ϋߓ_�p */
	Ek=Mk; for(i=0;i<10;i++) Ek=Mk+e*sin(Ek);	/* Kepler������ */

	/* ���Θ_�␳ */
	tr	=-2.0*sqrt(MUe)/C/C*e*sqrtA*sin(Ek);

	/* ���ԍ������߂� */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOC);
	tk	=tk0-psr/C;							/* �`�����Ԃ̕� */

	/* �q�����v�̕␳�ʂ��v�Z */
	dt	=get_ephemeris(prn,EPHM_AF0)
			+get_ephemeris(prn,EPHM_AF1)*tk
			+get_ephemeris(prn,EPHM_AF2)*tk*tk;

	return dt+tr-get_ephemeris(prn,EPHM_TGD);
}

/*------------------------------------------------------------
 * satellite_position() - �q���ʒu���v�Z
 *------------------------------------------------------------
 *  posxyz satellite_position(prn,wt,psr); �q���ʒu
 *    int prn;    �q��PRN�ԍ�(1�`)
 *    wtime wt;   �������w��
 *    double psr; �������̓`������[m]
 *------------------------------------------------------------
 *   ���O��set_ephemeris()�ɂ��G�t�F�����X���Z�b�g����Ă���
 * ����. ������get_ephemeris()���g�p����.
 *------------------------------------------------------------*/
posxyz satellite_position(int prn,wtime wt,double psr)
{
	int		i;
	double	tk,tk0,sqrtA,e,n,Ek,Mk,xk,yk,Omegak,
			vk,pk,uk,rk,ik,d_uk,d_rk,d_ik;
	posxyz	pos;

	/* ���ԍ������߂� */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOE);
	tk	=tk0-psr/C;							/* �`�����Ԃ̕� */

	/* ���S�ߓ_�pEk[rad]�����߂� */
	sqrtA	=get_ephemeris(prn,EPHM_sqrtA);
	e	=get_ephemeris(prn,EPHM_e); 		/* ���S�� */
	n	=sqrt(MUe)/sqrtA/sqrtA/sqrtA+get_ephemeris(prn,EPHM_d_n);
	Mk	=get_ephemeris(prn,EPHM_M0)+n*tk;	/* ���ϋߓ_�p */
	Ek=Mk; for(i=0;i<10;i++) Ek=Mk+e*sin(Ek);	/* Kepler������ */

	/* �O���ʓ��ɂ�����q���ʒu */
	rk	=sqrtA*sqrtA*(1.0-e*cos(Ek));		/* ���a�� */
	vk	=atan2((sqrt(1.0-e*e)*sin(Ek)),(cos(Ek)-e));/* �^�ߓ_�p */
	pk	=vk+get_ephemeris(prn,EPHM_omega);	/* �ܓx����[rad] */

	/* �␳�W����K�p���� */
	d_uk=get_ephemeris(prn,EPHM_Cus)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Cuc)*cos(2.0*pk);
	d_rk=get_ephemeris(prn,EPHM_Crs)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Crc)*cos(2.0*pk);
	d_ik=get_ephemeris(prn,EPHM_Cis)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Cic)*cos(2.0*pk);
	uk	=pk+d_uk;							/* �ܓx����[rad] */
	rk	=rk+d_rk;							/* ���a��[m] */
	ik	=get_ephemeris(prn,EPHM_i0)+d_ik
			+get_ephemeris(prn,EPHM_di)*tk;	/* �O���X�Ίp[rad] */

	/* �O���ʓ��ł̈ʒu */
	xk	=rk*cos(uk);
	yk	=rk*sin(uk);

	/* ����_�̌o�x[rad] */
	Omegak	=get_ephemeris(prn,EPHM_OMEGA0)
				+(get_ephemeris(prn,EPHM_dOmega)-dOMEGAe)*tk0
				-dOMEGAe*get_ephemeris(prn,EPHM_TOE);

	/* ECEF���W�n�ɕϊ� */
	pos.x	=xk*cos(Omegak)-yk*sin(Omegak)*cos(ik);
	pos.y	=xk*sin(Omegak)+yk*cos(Omegak)*cos(ik);
	pos.z	=yk*sin(ik);

	return pos;
}

/*------------------------------------------------------------
 * iono_correction() - �d���w�␳
 *------------------------------------------------------------
 *  double iono_correction(sat,usr,wt); �d���w�␳[m]
 *    posxyz sat; �q���ʒu
 *    posxyz usr; ���[�U�ʒu
 *    wtime wt;   �������w��
 *------------------------------------------------------------
 *   ���O��iono_alpha[],iono_beta[]���Z�b�g����Ă��邱��.
 * �Ԃ����l�͕␳�l�Ȃ̂ŁA�����ƂȂ邱�Ƃɒ���.
 *------------------------------------------------------------*/
#define	PER_MIN			(20.0*3600.0)
#define	NIGHT_DELAY		5.0e-9
#define	MAX_DELAY_TIME	(14.0*3600.0)
double iono_correction(posxyz sat,posxyz usr,wtime wt)
{
	int		i;
	double	x,y,f,lt,amp,per,phi_m,phi_i,lam_i,psi,el,az;
	posblh	usrblh;

	/* �p�ƕ��ʊp�����߂� */
	el	=rad_to_sc(elevation(sat,usr));
	az	=rad_to_sc(azimuth(sat,usr));

	/* �s�A�[�X�|�C���g�����߂� */
	usrblh	=xyz_to_blh(usr);
	psi		=0.0137/(el+0.11)-0.022;
	phi_i	=rad_to_sc(usrblh.lat)+psi*cos(az*PI);
	if (phi_i>0.416)	phi_i=0.416;
	if (phi_i<-0.416)	phi_i=-0.416;
	lam_i	=rad_to_sc(usrblh.lon)+psi*sin(az*PI)/cos(phi_i*PI);
	phi_m	=phi_i+0.064*cos((lam_i-1.617)*PI);

	/* �n���� */
	lt	=(double)SECONDS_DAY/2.0*lam_i+wt.sec;
	while(lt>SECONDS_DAY)	lt-=SECONDS_DAY;
	while(lt<0.0)			lt+=SECONDS_DAY;

	/* �R�T�C���֐��̎����ƐU�� */
	amp=0.0; per=0.0; for(i=0;i<IONO_PARAMETERS;i++) {
		amp	+=iono_alpha[i]*pow(phi_m,(double)i);
		per	+=iono_beta[i]*pow(phi_m,(double)i);
	}
	if (amp<0.0)		amp=0.0;
	if (per<PER_MIN)	per=PER_MIN;

	/* �X�ΌW�� */
	x	=0.53-el;
	f	=1.0+16.0*x*x*x;

	/* �␳�l�����߂ĕԂ� */
	x	=2.0*PI*(lt-MAX_DELAY_TIME)/per;
	while(x>PI)		x-=2.0*PI;
	while(x<-PI)	x+=2.0*PI;
	if (fabs(x)<1.57) {
		y	=amp*(1.0-x*x*(0.5-x*x/24.0));	/* ���� */
	} else y=0.0;							/* ��� */

	return -f*(NIGHT_DELAY+y)*C;		/* �����ɂ��ĕԂ� */
}

/*------------------------------------------------------------
 * tropo_correction() - �Η����␳
 *------------------------------------------------------------
 *  double tropo_correction(sat,usr); �Η����␳[m]
 *    posxyz sat; �q���ʒu
 *    posxyz usr; ���[�U�ʒu
 *------------------------------------------------------------
 *   �Ԃ����l�͕␳�l�Ȃ̂ŁA�����ƂȂ邱�Ƃɒ���.
 *------------------------------------------------------------*/
#define	TROPO_DELAY_ZENITH		2.47
#define	TROPO_SCALE_HEIGHT		(1.0/2.3E-5)
double tropo_correction(posxyz sat,posxyz usr)
{
	double	d;
	posblh	usrblh;

	usrblh	=xyz_to_blh(usr);
	if (usrblh.hgt<0.0) {
		d	=1.0;
	} else if (usrblh.hgt>TROPO_SCALE_HEIGHT) {
		d	=0.0;
	} else {
		d	=1.0-usrblh.hgt/TROPO_SCALE_HEIGHT;
	}

	return	-TROPO_DELAY_ZENITH*d*d*d*d*d
				/(sin(elevation(sat,usr))+0.0121);
}

/*------------------------------------------------------------
 * ���ʌv�Z
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * inverse_matrix() - �t�s����v�Z����
 *------------------------------------------------------------
 *  void inverse_matrix(a,n);
 *    double a[][]; ���̍s��(�t�s��ŏ㏑�������)
 *    int m;        �s��̎���(1�`MAX_M)
 *------------------------------------------------------------
 *   �^����ꂽ�s��̋t�s������߂�. ���ʂɂ��A���̍s��
 * �㏑�������.
 *------------------------------------------------------------*/
static void inverse_matrix(double a[MAX_M][MAX_M],int m)
{
	int		i,j,k;
	double	b[MAX_M][MAX_M+MAX_M];

	/* ����p�̍s������� */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			b[i][j]=a[i][j];
			if (i==j) b[i][j+m]=1.0; else b[i][j+m]=0.0;
		}
	}

	/* �K�E�X�̏����@ */
	for(i=0;i<m;i++) {
		/* ��i�s��b[i][i]�Ő��K������ */
		if (fabs(b[i][i])<=1E-10) {
			fprintf(stderr,"Cannot inverse matrix.\n");
			exit(2);
		}
		for(j=m+m-1;j>=i;j--) {
			b[i][j]/=b[i][i];
		}

		/* ���̍s�̑�i����������� */
		for(k=0;k<m;k++) if (k!=i) {
			for(j=m+m-1;j>=i;j--) {
				b[k][j]-=b[k][i]*b[i][j];
			}
		}
	}

	/* ���̍s����t�s��ŏ㏑������ */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			a[i][j]=b[i][j+m];
		}
	}
}

/*------------------------------------------------------------
 * compute_solution() - �ŏ����@�ŕ�����������
 *------------------------------------------------------------
 *  void compute_solution(G,dr,wgt,dx,cov,n,m);
 *    double G[][];   �f�U�C���s��(n�~m)
 *    double dr[];    �������̉E��(n��)
 *    double wgt[];   �d�݌W��(n��)/NULL:�d�݂Ȃ�
 *    double dx[];    �������̉��ŏ㏑�������(m��)
 *    double cov[][]; �����U�s��ŏ㏑�������(m�~m)
 *    int n;          �������̐�
 *    int m;          ���m���̐�
 *------------------------------------------------------------
 *   �^����ꂽ���������ŏ����@�ɂ�����. �d�݂��s�v�ȏꍇ
 * �� wgt=NULL �Ƃ��ČĂяo��.
 *------------------------------------------------------------*/
void compute_solution(double G[MAX_N][MAX_M],double dr[MAX_N],
	double wgt[MAX_N],double dx[MAX_M],double cov[MAX_M][MAX_M],
	int n,int m)
{
	int		i,j,k;
	double	w,a[MAX_M][MAX_N];

	/* GtG�����߂� */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			cov[i][j]=0.0;
			for(k=0;k<n;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[k];
				cov[i][j]+=G[k][i]*G[k][j]*w;
			}
		}
	}

	/* �t�s������߂�i���ꂪ�����U�s��C�ƂȂ�j */
	inverse_matrix(cov,m);

	/* Gt�������� */
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			a[i][j]=0.0;
			for(k=0;k<m;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[j];
				a[i][j]+=cov[i][k]*G[j][k]*w;
			}
		}
	}

	/* dr��������Ɖ��ɂȂ� */
	for(i=0;i<m;i++) {
		dx[i]=0.0;
		for(k=0;k<n;k++) {
			dx[i]+=a[i][k]*dr[k];
		}
	}
}

/*------------------------------------------------------------
 * _compute_position() - ���ʌv�Z�����̉��ʃ��[�`��
 *------------------------------------------------------------
 *  int _compute_position(wt,psr1,detail,sol,cov,dpsr,dpsr1,
 *      el,az); �q����
 *    wtime wt;       ��M�����i�^�C���X�^���v�j
 *    double psr1[];  L1�[������(PRN��; 0.0�Ȃ疳��)
 *    bool detail;    TRUE:�����v�Z
 *    double sol[];   �ʒu�E�N���b�N��(�����l�����ČĂяo��)
 *    double cov[][]; �����U�s��ŏ㏑�������(4�~4)
 *    double dpsr[];  �[�������̎c��(PRN��; 0.0�Ȃ疳��)
 *    double dpsr1[]; �c���̂����̑�C�x������(PRN��)
 *    double el[];    �q���p(PRN��; 0.0�Ȃ疳��)
 *    double az[];    �q�����ʊp(PRN��; 0.0�Ȃ疳��)
 *------------------------------------------------------------
 *   ���ʌv�Z�̌J�Ԃ��v�Z����񕪂������s����. psr1[]�ɂ�
 * PRN���ɋ[�����������ČĂяo������(�����ȉq���ɂ��Ă�
 * 0.0���Z�b�g���Ă���). sol[]�͈ʒu�ƃN���b�N�̏����l������
 * �Ăяo���ƁA�X�V���ꂽ�l���i�[�����(���e��������������).
 * cov[][]�͍ŏ����@�̋����U�s��Adpsr[]�͎c���x�N�g���A�܂�
 * el[],az[]�͉q���̋p�E���ʊp�ł��ꂼ��㏑�������.
 *------------------------------------------------------------*/
static int _compute_position(wtime wt,double *psr1,bool detail,
	double sol[MAX_M],double cov[MAX_M][MAX_M],double *dpsr,
	double *dpsr1,double *el,double *az)
{
	int		i,m,n=0,prn;
	double	r,satclk;
	double	G[MAX_N][MAX_M],dr[MAX_N],dx[MAX_M];
	double	wgt[MAX_N],cov2[MAX_M][MAX_M];
	posxyz	satpos,usrpos,pos,base;
	posenu	denu;
	posblh	blh;

	/* ��M�����i��M�@�N���b�N����␳����j */
	wt.sec	-=sol[3]/C;

	/* ���m���̐� */
	m=4;

	/* �b��̃��[�U�ʒu */
	usrpos.x=sol[0];
	usrpos.y=sol[1];
	usrpos.z=sol[2];

	/* �[���������L���ȉq���݂̂��g�p���� */
	for(prn=1;prn<=MAX_PRN;prn++) {
		dpsr[prn-1]	=0.0;
		dpsr1[prn-1]=0.0;
		el[prn-1]	=0.0;
		az[prn-1]	=0.0;
		if (psr1[prn-1]>0.0) {
			/* �q���N���b�N�Ɖq���ʒu���v�Z */
			r		=psr1[prn-1]-sol[3];
			satclk	=satellite_clock(prn,wt,r);
			r		=psr1[prn-1]-sol[3]+satclk*C;
			satclk	=satellite_clock(prn,wt,r);
			satpos	=satellite_position(prn,wt,r);
			el[prn-1]	=elevation(satpos,usrpos);
			az[prn-1]	=azimuth(satpos,usrpos);

			/* �f�U�C���s������� */
			r		=DIST(satpos,usrpos);
			G[n][0]	=-xyz_to_enu(satpos,usrpos).e/r;
			G[n][1]	=-xyz_to_enu(satpos,usrpos).n/r;
			G[n][2]	=-xyz_to_enu(satpos,usrpos).u/r;
			G[n][3]	=1.0;

			/* �[�������̏C���� */
			dr[n]	=psr1[prn-1]+satclk*C-(r+sol[3]);

			/* �d�݂�ݒ肷�� */
			if (el[prn-1]>0.0) {
				wgt[n]	=sin(el[prn-1])*sin(el[prn-1])/VAR_ZENITH;
			} else wgt[n]=0.0;
			if (wgt[n]<MIN_WEIGHT) wgt[n]=MIN_WEIGHT;

			/* �����v�Z */
			if (detail) {
				/* �p�}�X�N�ɂ��I�� */
				if (el[prn-1]<mask_angle) continue;

				/* �␳�l�������� */
				dpsr1[prn-1]+=iono_correction(satpos,usrpos,wt);
				dpsr1[prn-1]+=tropo_correction(satpos,usrpos);
				dr[n]		+=dpsr1[prn-1];
			}

			/* �c����Ԃ� */
			dpsr[prn-1]	=dr[n++];
		}
	}
	if (n<m) return 0;

	/* ������������ */
	compute_solution(G,dr,NULL,dx,cov,n,m);
	compute_solution(G,dr,wgt,dx,cov2,n,m);

	/* �����l�ɉ����� */
	denu.e	=dx[0];
	denu.n	=dx[1];
	denu.u	=dx[2];
	sol[0]	=enu_to_xyz(denu,usrpos).x;
	sol[1]	=enu_to_xyz(denu,usrpos).y;
	sol[2]	=enu_to_xyz(denu,usrpos).z;
	sol[3]	+=dx[m-1];

	return n;
}

/*------------------------------------------------------------
 * compute_position() - ���ʌv�Z����
 *------------------------------------------------------------
 *  void compute_position(wt,psr1);
 *    wtime wt;      ��M����
 *    double psr1[]; L1�[������
 *------------------------------------------------------------*/
#define	LOOP	8
void compute_position(wtime wt,double psr1[MAX_PRN])
{
	int		i,n,prn,loop,max_prn,latd,lond,latm,lonm,sum;
	bool	detail;
	double	sol[MAX_M],cov[MAX_M][MAX_M];
	double	dpsr[MAX_PRN],dpsr1[MAX_PRN],el[MAX_PRN],az[MAX_PRN];
	double	gdop,pdop,hdop,vdop,d,max_dpsr,latmf,lonmf;
	posxyz	pos;
	posblh	blh;
	posenu	enu;
	struct tm	tmbuf;

	/* �G�t�F�����X���Z�b�g */
	for(prn=1;prn<=MAX_PRN;prn++) {
		if (psr1[prn-1]>0.0) {
			if (!set_ephemeris(prn,wt,-1)) {
				/* �����ȉq�� */
				psr1[prn-1]=0.0;
			} else {
				/* health�t���O���`�F�b�N */
				if (get_ephemeris(prn,EPHM_health)!=0.0) {
					psr1[prn-1]=0.0;
				}
			}
		}
	}

	/* ���������� */
	for(i=0;i<MAX_M;i++) sol[i]=0.0;

	/* ���ʌv�Z */
	for(loop=0;loop<LOOP;loop++) {
		/* 3��ڂ܂ł͑e�v�Z�A����ȍ~�͐����v�Z */
		if (loop<3) detail=FALSE; else detail=TRUE;

		/* �v�Z���[�`�����Ăяo�� */
		if ((n=_compute_position(wt,psr1,detail,sol,cov,
			dpsr,dpsr1,el,az))<1) break;

		/* �c�����傫�ȉq���͏��� */
		if (loop>=LOOP-3) {
			max_prn	=0;
			for(prn=1;prn<=MAX_PRN;prn++) {
				if (psr1[prn-1]>0.0) {
					if (max_prn<1 || fabs(dpsr[prn-1])>max_dpsr) {
						max_prn	=prn;
						max_dpsr=fabs(dpsr[prn-1]);
					}
				}
			}
			if (max_prn>0 && max_dpsr>DPSR_THRESHOLD) {
				psr1[max_prn-1]=0.0;
			}
		}
	}
	if (n<1) return;

	/* ���ʂ��o�͂��� */
	pos.x	=sol[0];
	pos.y	=sol[1];
	pos.z	=sol[2];
	wt.sec	-=sol[3]/C;
	blh		=xyz_to_blh(pos);
	gdop	=sqrt(cov[0][0]+cov[1][1]+cov[2][2]+cov[3][3]);
	pdop	=sqrt(cov[0][0]+cov[1][1]+cov[2][2]);
	hdop	=sqrt(cov[0][0]+cov[1][1]);
	vdop	=sqrt(cov[2][2]);

	/* NMEA�t�H�[�}�b�g */
	if (output_nmea) {
		/* ��M���� */
		wt.sec	-=leap_sec;
		tmbuf	=wtime_to_date(wt);

		/* �o�ܓx��x�ƕ��ɕ����� */
		d		=fabs(rad_to_deg(blh.lat));
		latd	=(int)d;
		d		=(d-(double)latd)*60.0;
		latm	=(int)d/10;
		latmf	=d-(double)latm*10.0;
		d		=fabs(rad_to_deg(blh.lon));
		lond	=(int)d;
		d		=(d-(double)lond)*60.0;
		lonm	=(int)d/10;
		lonmf	=d-(double)lonm*10.0;

		/* �o�͓��e�𐶐����� */
		sprintf(linebuf,"$GPGGA,%2.2d%2.2d%2.2d,\
%d%d%.7f,%c,%d%d%.7f,%c,%d,%d,%.3f,%.3f,M,%.3f,M,,",
			tmbuf.tm_hour,			/* #2: ��M����[hhmmss] */
			tmbuf.tm_min,
			tmbuf.tm_sec,
			latd,latm,latmf,		/* #3: �ܓx[ddmm.mmm] */
			(blh.lat<0.0)?'S':'N',	/* #4: �ܓx�̕��� */
			lond,lonm,lonmf,		/* #5: �o�x[ddmm.mmm] */
			(blh.lon<0.0)?'W':'E',	/* #6: �o�x�̕��� */
			1,						/* #7: 1:�P��/2:�␳���� */
			n,						/* #8: �q���� */
			hdop,					/* #9: HDOP */
			blh.hgt-GEOIDAL_HEIGHT,	/* #10: �W��[m] */
			GEOIDAL_HEIGHT);		/* #12: �W�I�C�h��[m] */

		/* �`�F�b�N�T�� */
		sum=0; for(i=0;i<strlen(linebuf);i++) {
			sum^=linebuf[i];
		}

		/* �o�͂��� */
		printf("%s*%2.2X\n",linebuf,sum);

	/* ���Έʒu */
	} else if (output_rel) {
		/* ���Έʒu�ɕϊ� */
		enu=xyz_to_enu(pos,blh_to_xyz(rel_base));

		/* �o�͂��� */
		printf("%.3f,%.4f,%.4f,%.4f,%.4E,%d,%.3f,%.3f,%.3f,%.3f\n",
			wt.sec,				/* #1: ��M����[s] */
			enu.e,enu.n,enu.u,	/* #2-4: ��M�@�ʒu */
			sol[3]/C,			/* #5: ��M�@�N���b�N�덷[s] */
			n,					/* #6: �q���� */
			gdop,				/* #7: GDOP */
			pdop,				/* #8: PDOP */
			hdop,				/* #9: HDOP */
			vdop);				/* #10:VDOP */

	/* ECEF���W�ŏo�� */
	} else {
		printf("%.3f,%.4f,%.4f,%.4f,%.4E,%d,%.3f,%.3f,%.3f,%.3f\n",
			wt.sec,				/* #1: ��M����[s] */
			pos.x,pos.y,pos.z,	/* #2-4: ��M�@�ʒu */
			sol[3]/C,			/* #5: ��M�@�N���b�N�덷[s] */
			n,					/* #6: �q���� */
			gdop,				/* #7: GDOP */
			pdop,				/* #8: PDOP */
			hdop,				/* #9: HDOP */
			vdop);				/* #10:VDOP */
	}
}

/*------------------------------------------------------------
 * �ϑ��f�[�^�t�@�C���̏���
 *------------------------------------------------------------*/
/* RINEX�t�@�C���̏�� */
#define	MAX_OBSVS_LINE			5
#define	MAX_TYPES_OF_OBSV_LINE	9
#define	MAX_PRNS_LINE			12
#define	MAX_OBSVS				32

/* �ϑ��l�̎�� */
#define	MEAS_TYPES				4	/* ���̌��܂ł��L�� */
static char	*type_name[]={
	"C1","P2","L1","L2","P1",NULL
};	/* ���₷�ꍇ��P1�̑O�ɓ���AMEAS_TYPES�����₷ */
enum meas_type {
	TYPE_C1,TYPE_P2,TYPE_L1,TYPE_L2,TYPE_P1
};	/* type_name[]�ɑΉ������� */

/*------------------------------------------------------------
 * read_RINEX_OBS() - RINEX�ϑ��t�@�C����ǂݍ���
 *------------------------------------------------------------
 *  void read_RINEX_OBS(fp);
 *    FILE *fp; �ǂݍ��ރt�@�C��
 *------------------------------------------------------------*/
void read_RINEX_OBS(FILE *fp)
{
	int		i,j,n,prn,line,stat,num,lli,snr,obsvs=0,
			prn_list[MAX_PRN],type_num[MAX_OBSVS];
	bool	noerr=FALSE;
	double	d,obs[MEAS_TYPES][MAX_PRN];
	char	*p;
	wtime	wt;
	struct tm	tmbuf,tmbuf1={0,0,0,0,0,0};

	/* �w�b�_���� */
	fprintf(stderr,"Reading RINEX OBS... ");
	while(read_line(fp)) {
		if (is_comment("# / TYPES OF OBSERV")) {
			/* �ϑ��f�[�^�̐� */
			n=atoi(get_field(6));
			for(i=0;i<n;i++) {
				/* �K�v�Ȃ玟�̍s��ǂ� */
				if (i>0 && (i%MAX_TYPES_OF_OBSV_LINE)==0) {
					if (!read_line(fp)) goto ERROR;
					get_field(6);
				}
				/* �ϑ��f�[�^�� */
				p=get_field(6);
				while(isspace(*p)) p++;
				for(j=0;type_name[j]!=NULL;j++) {
					if (strcmp(p,type_name[j])==0) break;
				}
				type_num[i+obsvs]=j;	/* ���Ԃ��L������ */
			}
			obsvs+=n;
		} else if (is_comment("END OF HEADER")) break;
	}

	/* C1���Ȃ����P1�ő�p���� */
	for(i=0;i<obsvs;i++) {
		if (type_num[i]==TYPE_C1) break;
	}
	if (i>=obsvs) {
		for(i=0;i<obsvs;i++) {
			if (type_num[i]==TYPE_P1) type_num[i]=TYPE_C1;
		}
	}

	/* �{�̂�ǂ� */
	while(read_line(fp)) {
		/* �G�|�b�N */
		tmbuf.tm_year	=atoi(get_field(3));
		if (tmbuf.tm_year<80) tmbuf.tm_year+=100;
		tmbuf.tm_mon	=atoi(get_field(3))-1;
		tmbuf.tm_mday	=atoi(get_field(3));
		tmbuf.tm_hour	=atoi(get_field(3));
		tmbuf.tm_min	=atoi(get_field(3));
		d				=atof(get_field(11));
		tmbuf.tm_sec	=(int)(d+0.5);
		d				-=tmbuf.tm_sec;
		if (tmbuf.tm_mon>=0 && tmbuf.tm_mday>0) {
			wt		=date_to_wtime(tmbuf);
			wt.sec	+=d;
			/* �J�n���� */
			if (tmbuf1.tm_year<1) {
				fprintf(stderr,"%.2d/%.2d/%.2d %.2d:%.2d:%.2d - ",
					tmbuf.tm_year%100,tmbuf.tm_mon+1,tmbuf.tm_mday,
					tmbuf.tm_hour,tmbuf.tm_min,tmbuf.tm_sec);
			}
			tmbuf1	=tmbuf;
		} else wt.week=-1;

		/* �X�e�[�^�X�ƃ��R�[�h�� */
		stat		=atoi(get_field(3));
		num			=atoi(get_field(3));

		/* �C�x���g�t���O�̏��� */
		if (stat>=2) {
			/* �C�x���g�t���O���������|��\�� */
			fprintf(stderr,"Detected: EVENT FLAG %d\n",stat);
			if (stat>=6) continue;
			for(i=0;i<num*((obsvs-1)/MAX_OBSVS_LINE+1);i++) {
				if (!read_line(fp)) goto ERROR;
				get_field(60);
				if (isalpha(*get_field(1))) {
					i+=(obsvs-1)/MAX_OBSVS_LINE;
				}
			}
			continue;
		}
		/* �G�|�b�N��񂪂Ȃ��ꍇ */
		if (wt.week<0) {
			for(i=0;i<num*((obsvs-1)/MAX_OBSVS_LINE+1);i++) {
				if (!read_line(fp)) goto ERROR;
			}
			continue;
		}
		/* ���R�[�h�����`�F�b�N */
		if (num<1) continue;

		/* PRN�ԍ��̃��X�g��ǂݍ��� */
		for(i=0;i<num;i++) {
			/* �K�v�Ȃ玟�̍s��ǂ� */
			if (i>0 && (i%MAX_PRNS_LINE)==0) {
				if (!read_line(fp)) goto ERROR;
				get_field(32);
			}
			/* PRN�ԍ� */
			prn_list[i]=0;
			switch(*get_field(1)) {
				case ' ':
				case 'G':
					prn=atoi(get_field(2));
					if (prn>=1 && prn<=MAX_PRN) {
						prn_list[i]=prn;
					}
					break;
				default:
					get_field(2);
			}
		}

		/* ������ */
		for(i=0;i<MEAS_TYPES;i++) {
			for(prn=1;prn<=MAX_PRN;prn++) {
				obs[i][prn-1]=0.0;
			}
		}

		/* ���R�[�h��ǂݍ��� */
		for(line=0;line<num;line++) {
			if (!read_line(fp)) goto ERROR;
			for(i=0;i<obsvs;i++) {
				/* �K�v�Ȃ玟�̍s��ǂ� */
				if (i>0 && (i%MAX_OBSVS_LINE)==0) {
					if (!read_line(fp)) goto ERROR;
				}
				/* �ϑ��f�[�^ */
				prn	=prn_list[line];
				d	=atof(get_field(14));
				lli	=atoi(get_field(1));
				snr	=atoi(get_field(1));
				if (prn>0 && type_num[i]<MEAS_TYPES && d!=0.0) {
					obs[type_num[i]][prn-1]=d;
				}
			}
		}

		/* ���ʏ������Ăяo�� */
		compute_position(wt,obs[TYPE_C1]);
	}
	noerr=TRUE;

ERROR:
	/* �I������ */
	if (tmbuf1.tm_year>0) {
		fprintf(stderr,"%.2d/%.2d/%.2d %.2d:%.2d:%.2d\n",
			tmbuf1.tm_year%100,tmbuf1.tm_mon+1,tmbuf1.tm_mday,
			tmbuf1.tm_hour,tmbuf1.tm_min,tmbuf1.tm_sec);
	} else {
		fprintf(stderr,"No observation.\n");
	}

	/* �r���Ńt�@�C�����I����Ă����ꍇ */
	if (!noerr) {
		fprintf(stderr,"Error: Unexpected EOF.\n");
	}
}

/*------------------------------------------------------------
 * main() - ���C��
 *------------------------------------------------------------*/
void main(int argc,char **argv)
{
	int		i,opts;
	FILE	*fp;

	/* �I�v�V�����̏��� */
	opts=0;
	for(i=1;i<argc;i++) {
		if (argv[i][0]=='-') switch(argv[i][1]) {
			case 'n':
				output_nmea		=TRUE;
				opts++;
				break;
			case 'r':
				if (i+3<argc) {
					output_nmea	=FALSE;
					output_rel	=TRUE;
					rel_base.lat=deg_to_rad(atof(argv[++i]));
					rel_base.lon=deg_to_rad(atof(argv[++i]));
					rel_base.hgt=atof(argv[++i]);
					opts+=4;
				} else opts=argc;
				break;
			case 'm':
				if (i+1<argc) {
					mask_angle	=deg_to_rad(atof(argv[++i]));
					opts+=2;
				} else opts=argc;
				break;
			default:
				opts=argc;
		}
	}

	/* �R�}���h���C���������`�F�b�N */
	if (argc-opts<3) {
		printf("pos1 - standalone positioning.\n");
		printf("\n");
		printf("usage: pos1 [options] <nav> <obs>\n");
		printf("\n");
		printf("parameters:\n");
		printf("    <nav>       RINEX navigation file\n");
		printf("    <obs>       RINEX observation file\n");
		printf("\n");
		printf("options:\n");
		printf("    -n                   NMEA $GPGGA format\n");
		printf("    -r <lat> <lon> <hgt> Relative position\n");
		printf("    -m <mask>            Elevation mask [deg]\n");
		printf("\n");
		exit(0);
	}

	/* RINEX�t�@�C����ǂݍ��� */
	for(i=1+opts;i<argc;i++) {
		if ((fp=fopen(argv[i],"rt"))==NULL) {
			perror(argv[i]);
			exit(2);
		}
		if (toupper(argv[i][strlen(argv[i])-1])=='N') {
			/* RINEX�q�@�t�@�C�� */
			read_RINEX_NAV(fp);
		} else if (toupper(argv[i][strlen(argv[i])-1])=='O') {
			/* RINEX�ϑ��f�[�^�t�@�C�� */
			read_RINEX_OBS(fp);
		}
		fclose(fp);
	}

	exit(0);
}

