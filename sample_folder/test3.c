/*------------------------------------------------------------
 * TEST3.c - Practice for Position Computation.
 *------------------------------------------------------------*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*------------------------------------------------------------
 * 定数・構造体の定義
 *------------------------------------------------------------*/
/* 論理型 */
typedef	int	bool;
#define	TRUE	1
#define	FALSE	0

/* WGS-84定数 */
#define	PI		3.1415926535898		/* 円周率（IS-GPS-200） */
#define	C		2.99792458e8		/* 光速[m/s] */
#define	MUe		3.986005e14			/* 地球重力定数[m^3/s^2] */
#define	dOMEGAe	7.2921151467e-5		/* 地球自転角速度[rad/s] */
#define	Re		6378137.0			/* 地球半径[m] */
#define	Fe		(1.0/298.257223563)	/* 地球の扁平率 */

/* 角度の変換 */
#define	rad_to_deg(rad)		((rad)/PI*180.0)
#define	deg_to_rad(deg)		((deg)/180.0*PI)
#define	rad_to_sc(rad)		((rad)/PI)
#define	sc_to_rad(sc)		((sc)*PI)

/* 取り扱える行列の大きさ */
#define	MAX_N				16		/* 観測衛星数の上限 */
#define	MAX_M				4		/* 未知数の最大数 */
#define	MAX_PRN				32		/* 衛星番号の上限 */

/* 時間 */
#define	SECONDS_DAY			(3600L*24L)
#define	SECONDS_WEEK		(3600L*24L*7L)

/* 時刻を表す構造体 */
typedef struct {
	int		week;		/* 週番号 */
	double	sec;		/* 週初めからの経過時間[s] */
}	wtime;

/* 直交座標を表す構造体 */
typedef struct {
	double	x;			/* X座標[m] */
	double	y;			/* Y座標[m] */
	double	z;			/* Z座標[m] */
}	posxyz;
#define	SQ(x)			((x)*(x))
#define	DIST(a,b)		sqrt(SQ(a.x-b.x)+SQ(a.y-b.y)+SQ(a.z-b.z))

/* 経緯度を表す構造体 */
typedef struct {
	double	lat;		/* 緯度[rad] */
	double	lon;		/* 経度[rad] */
	double	hgt;		/* 高度（楕円体高）[m] */
}	posblh;

/* ENU座標を表す構造体 */
typedef struct {
	double	e;			/* East成分[m] */
	double	n;			/* North成分[m] */
	double	u;			/* Up成分[m] */
}	posenu;

/*------------------------------------------------------------
 * ファイルからの読込み
 *------------------------------------------------------------*/
/* 行バッファ */
#define	LINEBUF_LEN			256
static char	linebuf[LINEBUF_LEN],fieldbuf[LINEBUF_LEN];
static int	linepos			=0;

/* ファイルから行バッファに読み込む */
static bool read_line(FILE *fp)
{
	/* 行バッファをクリア */
	linepos		=0;
	linebuf[0]	='\0';

	/* ファイルから読み込む */
	return (fgets(linebuf,LINEBUF_LEN,fp)!=NULL);
}

/* 行バッファから文字列を取り出す（文字数指定、またはCSV形式） */
static char *get_field(int width)
{
	int		i;
	bool	quotef=TRUE;

	/* width==0ならCSV形式 */
	if (width<1) {
		width	=LINEBUF_LEN-1;
		quotef	=FALSE;
	}

	/* 指定された文字数あるいはコンマまでを読み取る */
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
	linepos		+=i;		/* 次の位置 */
	fieldbuf[i]	='\0';

	return fieldbuf;
}

/*------------------------------------------------------------
 * 座標変換
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * xyz_to_blh() - 直交座標から経緯度への変換
 *------------------------------------------------------------
 *  posblh xyz_to_blh(pos); 経緯度
 *    posxyz pos; 直交座標値
 *------------------------------------------------------------*/
posblh xyz_to_blh(posxyz pos)
{
	double	a,b,e,f,n,h,p,t,sint,cost;
	posblh	blh={0.0,0.0,-Re};

	/* 原点の場合 */
	if (pos.x==0.0 && pos.y==0.0 && pos.z==0.0) return blh;

	/* 楕円体のパラメータ */
	f	=Fe;				/* 扁平率 */
	a	=Re;				/* 長半径 */
	b	=a*(1.0-f);			/* 短半径 */
	e	=sqrt(f*(2.0-f));	/* 離心率 */

	/* 座標変換のためのパラメータ */
	h	=a*a-b*b;
	p	=sqrt(pos.x*pos.x+pos.y*pos.y);
	t	=atan2(pos.z*a,p*b);
	sint=sin(t);
	cost=cos(t);

	/* 経緯度への変換 */
	blh.lat =atan2(pos.z+h/b*sint*sint*sint,p-h/a*cost*cost*cost);
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	blh.lon	=atan2(pos.y,pos.x);
	blh.hgt	=(p/cos(blh.lat))-n;

	return blh;
}

/*------------------------------------------------------------
 * blh_to_xyz() - 経緯度から直交座標への変換
 *------------------------------------------------------------
 *  posxyz blh_to_xyz(blh); 直交座標値
 *    posblh blh; 経緯度
 *------------------------------------------------------------*/
posxyz blh_to_xyz(posblh blh)
{
	double	a,b,e,f,n;
	posxyz	pos;

	/* 楕円体のパラメータ */
	f	=Fe;				/* 扁平率 */
	a	=Re;				/* 長半径 */
	b	=a*(1.0-f);			/* 短半径 */
	e	=sqrt(f*(2.0-f));	/* 離心率 */

	/* 直交座標系への変換 */
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	pos.x	=(n+blh.hgt)*cos(blh.lat)*cos(blh.lon);
	pos.y	=(n+blh.hgt)*cos(blh.lat)*sin(blh.lon);
	pos.z	=(n*(1.0-e*e)+blh.hgt)*sin(blh.lat);

	return pos;
}

/*------------------------------------------------------------
 * xyz_to_enu() - 直交座標からENU座標への変換
 *------------------------------------------------------------
 *  posenu xyz_to_enu(pos,base); ENU座標値
 *    posxyz pos;  直交座標値
 *    posxyz base; 基準位置
 *------------------------------------------------------------*/
posenu xyz_to_enu(posxyz pos,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posenu	enu;

	/* 基準位置からの相対位置 */
	pos.x	-=base.x;
	pos.y	-=base.y;
	pos.z	-=base.z;

	/* 基準位置の経緯度 */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* 相対位置を回転させてENU座標に変換する */
	enu.e	=-pos.x*s1+pos.y*c1;
	enu.n	=-pos.x*c1*s2-pos.y*s1*s2+pos.z*c2;
	enu.u	=pos.x*c1*c2+pos.y*s1*c2+pos.z*s2;

	return enu;
}

/*------------------------------------------------------------
 * enu_to_xyz() - ENU座標から直交座標への変換
 *------------------------------------------------------------
 *  posxyz enu_to_xyz(enu,base); 直交座標値
 *    posenu enu;  ENU座標値
 *    posxyz base; 基準位置
 *------------------------------------------------------------*/
posxyz enu_to_xyz(posenu enu,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posxyz	pos;

	/* 基準位置の経緯度 */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* ENU座標を回転させて相対位置に変換する */
	pos.x	=-enu.e*s1-enu.n*c1*s2+enu.u*c1*c2;
	pos.y	=enu.e*c1-enu.n*s1*s2+enu.u*s1*c2;
	pos.z	=enu.n*c2+enu.u*s2;

	/* 基準位置に加える */
	pos.x	+=base.x;
	pos.y	+=base.y;
	pos.z	+=base.z;

	return pos;
}

/*------------------------------------------------------------
 * elevation() - 仰角を求める
 *------------------------------------------------------------
 *  double elevation(sat,usr); 仰角[rad]
 *    posxyz sat; 衛星の位置
 *    posxyz usr; ユーザ位置
 *------------------------------------------------------------*/
double elevation(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU座標に変換して仰角を求める */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.u,sqrt(enu.e*enu.e+enu.n*enu.n));
}

/*------------------------------------------------------------
 * azimuth() - 方位角を求める
 *------------------------------------------------------------
 *  double azimuth(sat,usr); 方位角[rad]
 *    posxyz sat; 衛星の位置
 *    posxyz usr; ユーザ位置
 *------------------------------------------------------------*/
double azimuth(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU座標に変換して方位角を求める */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.e,enu.n);
}

/*------------------------------------------------------------
 * 時刻の変換
 *------------------------------------------------------------*/
/* カレンダ値の開始年 */
#define	TIME_T_BASE_YEAR	1970

/* 1980年1月6日00:00:00のカレンダ値 */
#define	TIME_T_ORIGIN		315964800L

/* mktime()関数のGMT版（gmtime()関数に対応） */
static time_t mktime2(struct tm *tm)
{
	int		i;
	long	days=0L;
	static int	days_month[]={
		31,28,31,30,31,30,31,31,30,31,30,31
	};

	/* 経過日数を得る */
	for(i=TIME_T_BASE_YEAR;i<tm->tm_year+1900;i++) {
		days+=(i%4==0)?366:365;
	}
	for(i=1;i<tm->tm_mon+1;i++) {
		days+=days_month[i-1];
		if (i==2 && tm->tm_year%4==0) days++;
	}
	days+=tm->tm_mday-1;

	/* カレンダ値を返す */
	return ((days*24+tm->tm_hour)*60+tm->tm_min)*60+tm->tm_sec;
}

/*------------------------------------------------------------
 * wtime_to_date() - 週番号・秒から日時への変換
 *------------------------------------------------------------
 *  struct tm wtime_to_date(wt); 日時への変換結果
 *    wtime wt; 週番号・秒
 *------------------------------------------------------------*/
struct tm wtime_to_date(wtime wt)
{
	time_t	t;

	/* 基準日からの経過時間を加えて，カレンダ値を得る */
	t	=(long)wt.week*SECONDS_WEEK+TIME_T_ORIGIN
			+(long)((wt.sec>0.0)?wt.sec+0.5:wt.sec-0.5);

	return *gmtime(&t);
}

/*------------------------------------------------------------
 * date_to_wtime() - 日時から週番号・秒への変換
 *------------------------------------------------------------
 *  wtime date_to_wtime(tmbuf); 週番号・秒への変換結果
 *    struct tm tmbuf; 日時を指定
 *------------------------------------------------------------*/
wtime date_to_wtime(struct tm tmbuf)
{
	time_t	t;
	wtime	wt;

	/* 指定された時刻のカレンダ値 */
	t	=mktime2(&tmbuf);

	/* 1週間の秒数で割った商と余り */
	wt.week	=(t-TIME_T_ORIGIN)/SECONDS_WEEK;
	wt.sec	=(t-TIME_T_ORIGIN)%SECONDS_WEEK;

	return wt;
}

/*------------------------------------------------------------
 * 航法メッセージの処理
 *------------------------------------------------------------*/
/* RINEXファイルの情報 */
#define	RINEX_POS_COMMENT		60
#define	RINEX_NAV_LINES			8
#define	RINEX_NAV_FIELDS_LINE	4

/* 記憶するエフェメリスの最大数 */
#define	MAX_EPHMS				20

/* エフェメリスの有効期限[h] */
#define	EPHEMERIS_EXPIRE		2.0

/* エフェメリスを格納するための構造体 */
typedef	struct {
	int		week;		/* 週番号 */
	double	data[RINEX_NAV_LINES*RINEX_NAV_FIELDS_LINE];
}	ephm_info;

/* パラメータ番号を定義 */
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

/* エフェメリスを保持する配列 */
static ephm_info	ephm_buf[MAX_PRN][MAX_EPHMS];
static int			ephm_count[MAX_PRN];
static int			current_ephm[MAX_PRN];
static int			current_week=-1;

/* 閏秒の情報 */
static int			leap_sec	=0;

/* 電離層補正情報 */
#define	IONO_PARAMETERS			4
static double	iono_alpha[IONO_PARAMETERS];
static double	iono_beta[IONO_PARAMETERS];

/* 文字列→実数の変換 */
static double atof2(char *str)
{
	char	*p;

	/* 'D'を'E'に変換する */
	for(p=str;*p!='\0';p++) {
		if (*p=='D' || *p=='d') *p='E';
	}

	/* 実数に変換して返す */
	return atof(str);
}

/* コメント情報を調べる */
static bool is_comment(char *str)
{
	return (strncmp(linebuf+RINEX_POS_COMMENT,str,strlen(str))==0);
}

/*------------------------------------------------------------
 * read_RINEX_NAV() - RINEX航法ファイルを読み込む
 *------------------------------------------------------------
 *  void read_RINEX_NAV(fp);
 *    FILE *fp; 読み込むファイル
 *------------------------------------------------------------*/
void read_RINEX_NAV(FILE *fp)
{
	int		i,j,n,prn,line;
	bool	noerr=FALSE;
	double	d,t;
	wtime	wt;
	struct tm	tmbuf;
	ephm_info	info;

	/* 初期化 */
	if (current_week<0) {
		for(i=0;i<IONO_PARAMETERS;i++) {
			iono_alpha[i]	=0.0;
			iono_beta[i]	=0.0;
		}
		for(i=0;i<MAX_PRN;i++) {
			ephm_count[i]	=0;
		}
	}

	/* ヘッダ部分 */
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

	/* 本体を読む */
	fprintf(stderr,"Reading RINEX NAV... ");
	while(read_line(fp)) {
		n=0;
		for(line=0;line<RINEX_NAV_LINES;line++) {
			/* 最初のデータを読み込む */
			if (line==0) {
				/* 最初の行 */
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
				/* 2行目以降 */
				if (!read_line(fp)) goto ERROR;
				d	=atof2(get_field(22));
			}
			info.data[n++]=d;

			/* 残りのデータ */
			for(i=1;i<RINEX_NAV_FIELDS_LINE;i++) {
				d	=atof2(get_field(19));
				info.data[n++]=d;
			}
		}
		if (prn<1) continue;

		/* 週番号を揃える */
		t	=(info.week-info.data[EPHM_WEEK])*SECONDS_WEEK;
		info.week			=info.data[EPHM_WEEK];
		info.data[EPHM_TOC]	+=t;
		current_week		=info.week;

		/* すでに同じものがないかどうか */
		for(i=0;i<ephm_count[prn-1];i++) {
			/* 週番号が一致するものが対象 */
			if (ephm_buf[prn-1][i].week!=info.week) continue;

			/* IODCが一致しているか */
			if (ephm_buf[prn-1][i].data[EPHM_IODC]
					==info.data[EPHM_IODC]) {
				/* 送信時刻の早いものを残す */
				if (info.data[EPHM_TOT]
						<ephm_buf[prn-1][i].data[EPHM_TOT]) {
					ephm_buf[prn-1][i]=info;
				}
				prn=0;
				break;
			}
		}
		if (prn<1) continue;

		/* 配列に格納する（送信時刻の昇順） */
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
	/* エフェメリスがない場合 */
	if (current_week<0) {
		fprintf(stderr,"Error: No ephemeris information.\n");
		return;
	}

	/* エフェメリスのある衛星の数 */
	n=0; for(prn=1;prn<=MAX_PRN;prn++) {
		if (ephm_count[prn-1]>0) n++;
	}
	fprintf(stderr,"week %d: %d satellites\n",current_week,n);

	/* 途中でファイルが終わっていた場合 */
	if (!noerr) {
		fprintf(stderr,"Error: Unexpected EOF.\n");
	}
}

/*------------------------------------------------------------
 * set_ephemeris() - エフェメリスをセット
 *------------------------------------------------------------
 *  bool set_ephemeris(prn,wt,iode); TRUE:セットした
 *    int prn;  衛星PRN番号(1～)
 *    wtime wt; 時刻を指定
 *    int iode; IODEを指定/-1:指定なし
 *------------------------------------------------------------*/
bool set_ephemeris(int prn,wtime wt,int iode)
{
	int		i;
	double	t,t0;

	/* 新しいエフェメリスから探す */
	for(i=ephm_count[prn-1]-1;i>=0;i--) {
		/* 有効期限内であること */
		t0	=(ephm_buf[prn-1][i].week-wt.week)*SECONDS_WEEK;
		t	=ephm_buf[prn-1][i].data[EPHM_TOC]+t0;
		if (wt.sec<t-EPHEMERIS_EXPIRE*3600.0-0.1 ||
				t+EPHEMERIS_EXPIRE*3600.0+0.1<wt.sec) continue;

		/* IODEをチェック */
		if (iode>=0) {
			if (ephm_buf[prn-1][i].data[EPHM_IODE]==iode) break;
			continue;
		}

		/* 受信時刻 */
		t	=ephm_buf[prn-1][i].data[EPHM_TOT]+t0;
		if (t<wt.sec+0.1) break;	/* wtより前の時刻 */
	}

	/* カレント情報としてセットする */
	current_ephm[prn-1]=i;
	return (i>=0);
}

/*------------------------------------------------------------
 * get_ephemeris() - エフェメリスのパラメータを得る
 *------------------------------------------------------------
 *  double get_ephemeris(prn,para); パラメータ値
 *    int prn;  衛星PRN番号(1～)
 *    int para; パラメータ番号(0～)
 *------------------------------------------------------------
 *   事前にset_ephemeris()によりエフェメリスがセットされている
 * こと.
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
 * satellite_clock() - 衛星クロック誤差を計算
 *------------------------------------------------------------
 *  double satellite_clock(prn,wt,psr); 衛星クロック誤差[s]
 *    int prn;    衛星PRN番号(1～)
 *    wtime wt;   時刻を指定
 *    double psr; 見かけの伝搬距離[m]
 *------------------------------------------------------------
 *   事前にset_ephemeris()によりエフェメリスがセットされている
 * こと. 内部でget_ephemeris()を使用する.
 *------------------------------------------------------------*/
double satellite_clock(int prn,wtime wt,double psr)
{
	int		i;
	double	tk,tk0,sqrtA,e,n,Ek,Mk,dt,tr;

	/* 時間差を求める */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOE);
	tk	=tk0-psr/C;							/* 伝搬時間の分 */

	/* 離心近点角Ek[rad]を求める */
	sqrtA	=get_ephemeris(prn,EPHM_sqrtA);
	e	=get_ephemeris(prn,EPHM_e); 		/* 離心率 */
	n	=sqrt(MUe)/sqrtA/sqrtA/sqrtA+get_ephemeris(prn,EPHM_d_n);
	Mk	=get_ephemeris(prn,EPHM_M0)+n*tk;	/* 平均近点角 */
	Ek=Mk; for(i=0;i<10;i++) Ek=Mk+e*sin(Ek);	/* Kepler方程式 */

	/* 相対論補正 */
	tr	=-2.0*sqrt(MUe)/C/C*e*sqrtA*sin(Ek);

	/* 時間差を求める */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOC);
	tk	=tk0-psr/C;							/* 伝搬時間の分 */

	/* 衛星時計の補正量を計算 */
	dt	=get_ephemeris(prn,EPHM_AF0)
			+get_ephemeris(prn,EPHM_AF1)*tk
			+get_ephemeris(prn,EPHM_AF2)*tk*tk;

	return dt+tr-get_ephemeris(prn,EPHM_TGD);
}

/*------------------------------------------------------------
 * satellite_position() - 衛星位置を計算
 *------------------------------------------------------------
 *  posxyz satellite_position(prn,wt,psr); 衛星位置
 *    int prn;    衛星PRN番号(1～)
 *    wtime wt;   時刻を指定
 *    double psr; 見かけの伝搬距離[m]
 *------------------------------------------------------------
 *   事前にset_ephemeris()によりエフェメリスがセットされている
 * こと. 内部でget_ephemeris()を使用する.
 *------------------------------------------------------------*/
posxyz satellite_position(int prn,wtime wt,double psr)
{
	int		i;
	double	tk,tk0,sqrtA,e,n,Ek,Mk,xk,yk,Omegak,
			vk,pk,uk,rk,ik,d_uk,d_rk,d_ik;
	posxyz	pos;

	/* 時間差を求める */
	tk0	=(wt.week-get_ephemeris(prn,EPHM_WEEK))*SECONDS_WEEK
			+wt.sec-get_ephemeris(prn,EPHM_TOE);
	tk	=tk0-psr/C;							/* 伝搬時間の分 */

	/* 離心近点角Ek[rad]を求める */
	sqrtA	=get_ephemeris(prn,EPHM_sqrtA);
	e	=get_ephemeris(prn,EPHM_e); 		/* 離心率 */
	n	=sqrt(MUe)/sqrtA/sqrtA/sqrtA+get_ephemeris(prn,EPHM_d_n);
	Mk	=get_ephemeris(prn,EPHM_M0)+n*tk;	/* 平均近点角 */
	Ek=Mk; for(i=0;i<10;i++) Ek=Mk+e*sin(Ek);	/* Kepler方程式 */

	/* 軌道面内における衛星位置 */
	rk	=sqrtA*sqrtA*(1.0-e*cos(Ek));		/* 動径長 */
	vk	=atan2((sqrt(1.0-e*e)*sin(Ek)),(cos(Ek)-e));/* 真近点角 */
	pk	=vk+get_ephemeris(prn,EPHM_omega);	/* 緯度引数[rad] */

	/* 補正係数を適用する */
	d_uk=get_ephemeris(prn,EPHM_Cus)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Cuc)*cos(2.0*pk);
	d_rk=get_ephemeris(prn,EPHM_Crs)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Crc)*cos(2.0*pk);
	d_ik=get_ephemeris(prn,EPHM_Cis)*sin(2.0*pk)
			+get_ephemeris(prn,EPHM_Cic)*cos(2.0*pk);
	uk	=pk+d_uk;							/* 緯度引数[rad] */
	rk	=rk+d_rk;							/* 動径長[m] */
	ik	=get_ephemeris(prn,EPHM_i0)+d_ik
			+get_ephemeris(prn,EPHM_di)*tk;	/* 軌道傾斜角[rad] */

	/* 軌道面内での位置 */
	xk	=rk*cos(uk);
	yk	=rk*sin(uk);

	/* 昇交点の経度[rad] */
	Omegak	=get_ephemeris(prn,EPHM_OMEGA0)
				+(get_ephemeris(prn,EPHM_dOmega)-dOMEGAe)*tk0
				-dOMEGAe*get_ephemeris(prn,EPHM_TOE);

	/* ECEF座標系に変換 */
	pos.x	=xk*cos(Omegak)-yk*sin(Omegak)*cos(ik);
	pos.y	=xk*sin(Omegak)+yk*cos(Omegak)*cos(ik);
	pos.z	=yk*sin(ik);

	return pos;
}

/*------------------------------------------------------------
 * iono_correction() - 電離層補正
 *------------------------------------------------------------
 *  double iono_correction(sat,usr,wt); 電離層補正[m]
 *    posxyz sat; 衛星位置
 *    posxyz usr; ユーザ位置
 *    wtime wt;   時刻を指定
 *------------------------------------------------------------
 *   事前にiono_alpha[],iono_beta[]がセットされていること.
 * 返される値は補正値なので、負数となることに注意.
 *------------------------------------------------------------*/
#define	PER_MIN			(20.0*3600.0)
#define	NIGHT_DELAY		5.0e-9
#define	MAX_DELAY_TIME	(14.0*3600.0)
double iono_correction(posxyz sat,posxyz usr,wtime wt)
{
	int		i;
	double	x,y,f,lt,amp,per,phi_m,phi_i,lam_i,psi,el,az;
	posblh	usrblh;

	/* 仰角と方位角を求める */
	el	=rad_to_sc(elevation(sat,usr));
	az	=rad_to_sc(azimuth(sat,usr));

	/* ピアースポイントを求める */
	usrblh	=xyz_to_blh(usr);
	psi		=0.0137/(el+0.11)-0.022;
	phi_i	=rad_to_sc(usrblh.lat)+psi*cos(az*PI);
	if (phi_i>0.416)	phi_i=0.416;
	if (phi_i<-0.416)	phi_i=-0.416;
	lam_i	=rad_to_sc(usrblh.lon)+psi*sin(az*PI)/cos(phi_i*PI);
	phi_m	=phi_i+0.064*cos((lam_i-1.617)*PI);

	/* 地方時 */
	lt	=(double)SECONDS_DAY/2.0*lam_i+wt.sec;
	while(lt>SECONDS_DAY)	lt-=SECONDS_DAY;
	while(lt<0.0)			lt+=SECONDS_DAY;

	/* コサイン関数の周期と振幅 */
	amp=0.0; per=0.0; for(i=0;i<IONO_PARAMETERS;i++) {
		amp	+=iono_alpha[i]*pow(phi_m,(double)i);
		per	+=iono_beta[i]*pow(phi_m,(double)i);
	}
	if (amp<0.0)		amp=0.0;
	if (per<PER_MIN)	per=PER_MIN;

	/* 傾斜係数 */
	x	=0.53-el;
	f	=1.0+16.0*x*x*x;

	/* 補正値を求めて返す */
	x	=2.0*PI*(lt-MAX_DELAY_TIME)/per;
	while(x>PI)		x-=2.0*PI;
	while(x<-PI)	x+=2.0*PI;
	if (fabs(x)<1.57) {
		y	=amp*(1.0-x*x*(0.5-x*x/24.0));	/* 昼間 */
	} else y=0.0;							/* 夜間 */

	return -f*(NIGHT_DELAY+y)*C;		/* 距離にして返す */
}

/*------------------------------------------------------------
 * tropo_correction() - 対流圏補正
 *------------------------------------------------------------
 *  double tropo_correction(sat,usr); 対流圏補正[m]
 *    posxyz sat; 衛星位置
 *    posxyz usr; ユーザ位置
 *------------------------------------------------------------
 *   返される値は補正値なので、負数となることに注意.
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
 * 測位計算
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * inverse_matrix() - 逆行列を計算する
 *------------------------------------------------------------
 *  void inverse_matrix(a,n);
 *    double a[][]; 元の行列(逆行列で上書きされる)
 *    int m;        行列の次元(1～MAX_M)
 *------------------------------------------------------------
 *   与えられた行列の逆行列を求める. 結果により、元の行列が
 * 上書きされる.
 *------------------------------------------------------------*/
static void inverse_matrix(double a[MAX_M][MAX_M],int m)
{
	int		i,j,k;
	double	b[MAX_M][MAX_M+MAX_M];

	/* 操作用の行列をつくる */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			b[i][j]=a[i][j];
			if (i==j) b[i][j+m]=1.0; else b[i][j+m]=0.0;
		}
	}

	/* ガウスの消去法 */
	for(i=0;i<m;i++) {
		/* 第i行をb[i][i]で正規化する */
		if (fabs(b[i][i])<=1E-10) {
			fprintf(stderr,"Cannot inverse matrix.\n");
			exit(2);
		}
		for(j=m+m-1;j>=i;j--) {
			b[i][j]/=b[i][i];
		}

		/* 他の行の第i列を消去する */
		for(k=0;k<m;k++) if (k!=i) {
			for(j=m+m-1;j>=i;j--) {
				b[k][j]-=b[k][i]*b[i][j];
			}
		}
	}

	/* 元の行列を逆行列で上書きする */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			a[i][j]=b[i][j+m];
		}
	}
}

/*------------------------------------------------------------
 * compute_solution() - 最小二乗法で方程式を解く
 *------------------------------------------------------------
 *  void compute_solution(G,dr,wgt,dx,cov,n,m);
 *    double G[][];   デザイン行列(n×m)
 *    double dr[];    方程式の右辺(n次)
 *    double wgt[];   重み係数(n次)/NULL:重みなし
 *    double dx[];    方程式の解で上書きされる(m次)
 *    double cov[][]; 共分散行列で上書きされる(m×m)
 *    int n;          方程式の数
 *    int m;          未知数の数
 *------------------------------------------------------------
 *   与えられた方程式を最小二乗法により解く. 重みが不要な場合
 * は wgt=NULL として呼び出す.
 *------------------------------------------------------------*/
void compute_solution(double G[MAX_N][MAX_M],double dr[MAX_N],
	double wgt[MAX_N],double dx[MAX_M],double cov[MAX_M][MAX_M],
	int n,int m)
{
	int		i,j,k;
	double	w,a[MAX_M][MAX_N];

	/* GtGを求める */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			cov[i][j]=0.0;
			for(k=0;k<n;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[k];
				cov[i][j]+=G[k][i]*G[k][j]*w;
			}
		}
	}

	/* 逆行列を求める（これが共分散行列Cとなる） */
	inverse_matrix(cov,m);

	/* Gtをかける */
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			a[i][j]=0.0;
			for(k=0;k<m;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[j];
				a[i][j]+=cov[i][k]*G[j][k]*w;
			}
		}
	}

	/* drをかけると解になる */
	for(i=0;i<m;i++) {
		dx[i]=0.0;
		for(k=0;k<n;k++) {
			dx[i]+=a[i][k]*dr[k];
		}
	}
}

/*------------------------------------------------------------
 * _compute_position() - 測位計算処理の下位ルーチン
 *------------------------------------------------------------
 *  int _compute_position(wt,psr1,detail,sol,cov,dpsr,dpsr1,
 *      el,az); 衛星数
 *    wtime wt;       受信時刻（タイムスタンプ）
 *    double psr1[];  L1擬似距離(PRN順; 0.0なら無効)
 *    bool detail;    TRUE:精密計算
 *    double sol[];   位置・クロック解(初期値を入れて呼び出す)
 *    double cov[][]; 共分散行列で上書きされる(4×4)
 *    double dpsr[];  擬似距離の残差(PRN順; 0.0なら無効)
 *    double dpsr1[]; 残差のうちの大気遅延成分(PRN順)
 *    double el[];    衛星仰角(PRN順; 0.0なら無効)
 *    double az[];    衛星方位角(PRN順; 0.0なら無効)
 *------------------------------------------------------------
 *   測位計算の繰返し計算を一回分だけ実行する. psr1[]には
 * PRN順に擬似距離を入れて呼び出すこと(無効な衛星については
 * 0.0をセットしておく). sol[]は位置とクロックの初期値を入れて
 * 呼び出すと、更新された値が格納される(内容が書き換えられる).
 * cov[][]は最小二乗法の共分散行列、dpsr[]は残差ベクトル、また
 * el[],az[]は衛星の仰角・方位角でそれぞれ上書きされる.
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

	/* 受信時刻（受信機クロック分を補正する） */
	wt.sec	-=sol[3]/C;

	/* 未知数の数 */
	m=4;

	/* 暫定のユーザ位置 */
	usrpos.x=sol[0];
	usrpos.y=sol[1];
	usrpos.z=sol[2];

	/* 擬似距離が有効な衛星のみを使用する */
	for(prn=1;prn<=MAX_PRN;prn++) {
		dpsr[prn-1]	=0.0;
		dpsr1[prn-1]=0.0;
		el[prn-1]	=0.0;
		az[prn-1]	=0.0;
		if (psr1[prn-1]>0.0) {
			/* 衛星クロックと衛星位置を計算 */
			r		=psr1[prn-1]-sol[3];
			satclk	=satellite_clock(prn,wt,r);
			r		=psr1[prn-1]-sol[3]+satclk*C;
			satclk	=satellite_clock(prn,wt,r);
			satpos	=satellite_position(prn,wt,r);
			el[prn-1]	=elevation(satpos,usrpos);
			az[prn-1]	=azimuth(satpos,usrpos);

			/* デザイン行列をつくる */
			r		=DIST(satpos,usrpos);
			G[n][0]	=(sol[0]-satpos.x)/r;
			G[n][1]	=(sol[1]-satpos.y)/r;
			G[n][2]	=(sol[2]-satpos.z)/r;
			G[n][3]	=1.0;

			/* 擬似距離の修正量 */
			dr[n]	=psr1[prn-1]+satclk*C-(r+sol[3]);

			/* 精密計算 */
			if (detail) {
				/* 補正値を加える */
				dpsr1[prn-1]+=iono_correction(satpos,usrpos,wt);
				dpsr1[prn-1]+=tropo_correction(satpos,usrpos);
				dr[n]		+=dpsr1[prn-1];
			}

			/* 残差を返す */
			dpsr[prn-1]	=dr[n++];
		}
	}
	if (n<m) return 0;

	/* 方程式を解く */
	compute_solution(G,dr,NULL,dx,cov,n,m);
	/* 後で使うために空けてある */

	/* 初期値に加える */
	for(i=0;i<4;i++) {
		sol[i]+=dx[i];
	}
	/* 後で使うために空けてある
	 *
	 *
	 * */

	return n;
}

/*------------------------------------------------------------
 * compute_position() - 測位計算処理
 *------------------------------------------------------------
 *  void compute_position(wt,psr1);
 *    wtime wt;      受信時刻
 *    double psr1[]; L1擬似距離
 *------------------------------------------------------------*/
#define	LOOP	8
void compute_position(wtime wt,double psr1[MAX_PRN])
{
	int		i,n,prn,loop,max_prn,latd,lond,latm,lonm,sum;
	bool	detail;
	double	sol[MAX_M],cov[MAX_M][MAX_M];
	double	dpsr[MAX_PRN],dpsr1[MAX_PRN],el[MAX_PRN],az[MAX_PRN];

	/* エフェメリスをセット */
	for(prn=1;prn<=MAX_PRN;prn++) {
		if (psr1[prn-1]>0.0) {
			if (!set_ephemeris(prn,wt,-1)) {
				/* 無効な衛星 */
				psr1[prn-1]=0.0;
			} else {
				/* healthフラグをチェック */
				if (get_ephemeris(prn,EPHM_health)!=0.0) {
					psr1[prn-1]=0.0;
				}
			}
		}
	}

	/* 解を初期化 */
	for(i=0;i<MAX_M;i++) sol[i]=0.0;

	/* 測位計算 */
	for(loop=0;loop<LOOP;loop++) {
		/* 3回目までは粗計算、それ以降は精密計算 */
		if (loop<3) detail=FALSE; else detail=TRUE;

		/* 計算ルーチンを呼び出す */
		if ((n=_compute_position(wt,psr1,detail,sol,cov,
			dpsr,dpsr1,el,az))<1) break;

		/* 途中経過を出力する */
		printf("LOOP %d: X=%.3f, Y=%.3f, Z=%.3f, s=%.4E\n",
			loop+1,sol[0],sol[1],sol[2],sol[3]/C);
	}
}

/*------------------------------------------------------------
 * main() - メイン
 *------------------------------------------------------------*/
#define	SATS	5
static int		prn[SATS]={
	5,14,16,22,25,
};
static double	range[SATS]={
	23545777.534,	/* PRN 05 */
	20299789.570,	/* PRN 14 */
	24027782.537,	/* PRN 16 */
	24367716.061,	/* PRN 22 */
	22169926.127,	/* PRN 25 */
};

void main(int argc,char **argv)
{
	int		i;
	double	psr1[MAX_PRN];
	wtime	wt;
	FILE	*fp;

	/* RINEX航法ファイルを読み込む */
	if (argc<2) {
		fprintf(stderr,"test3 <RINEX-NAV>\n");
		exit(0);
	} else if ((fp=fopen(argv[1],"rt"))==NULL) {
		perror(argv[1]);
		exit(2);
	} else {
		read_RINEX_NAV(fp);
		fclose(fp);
	}

	/* 時刻を指定 */
	wt.week	=1349;		/* 05/11/13～19の週 */
	wt.sec	=86400.0;	/* 月曜日の00:00:00 */

	/* 擬似距離をセット */
	for(i=0;i<MAX_PRN;i++) {
		psr1[i]=0.0;
	}
	for(i=0;i<SATS;i++) {
		psr1[prn[i]-1]=range[i];
	}

	/* 測位計算を実行する */
	compute_position(wt,psr1);

	exit(0);
}

