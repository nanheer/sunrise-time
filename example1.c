#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eph_manager.h"
#include "novas.h"
#include "novascon.h"
#include<time.h>

int main (void)
{
    
    // UTC date and time
	time_t timep;
	struct tm *p;
	time(&timep);
	p = gmtime(&timep);
	// UTC date and time
	short int year = 1900 + p->tm_year; //UTC
	short int month = 1 + p->tm_mon;
	short int day = p->tm_mday;
	double hour = p->tm_hour + p->tm_min / 60.0 + p->tm_sec / 3600.0;//hour,格林尼治时间而不是北京时间
    const short int leap_secs = 37; //sec 1921=0.0
	int a, b;
	double c;
	int flag = 0;//flag=0计算默认计算机时间
	if (flag == 1){//flag=0计算手动输入时间
		year = 2018; //UTC
		month = 11;
		day = 20;
		hour =20.0;
	}
    //EOP
    const double ut1_utc = 0.06723; /* UT1-UTC (s)  */
    const double x_pole = 0.1538; /* polar motion (arcsec) */
    const double y_pole = 0.4389;

    //accuracy
   const short int accuracy = 0; //0 ... full accuracy, 1 ... reduced accuracy equ2hor
   short int error = 0;
   short int de_num = 0;
    short int coord_sys=1;  //0 ... GCRS or "local GCRS", 1 ...True equator and equinox of date, = 3 ... astrometric coordinates, i.e., without light deflection or aberration.
    
    double ddzz=0.00084; //deg 0.00084
    double ddsec=0.001; //sec 0.001

	//printf("北京时间or格林尼治时间：%f", hour);
    // observer on the surface of the Earth
	const double latitude = 31.0913;;  /* latitude (degree) */
   const double longitude = 121.16; /* east longitude (degree) */
	//const double latitude = 34.17625;;  /* latitude (degree) */
	//const double longitude = 115.928; /* east longitude (degree) */
	const double height = 12.2; /* height a.s.l (m) */
   const double temperature = 25.0; /* temperature T (K=273.15+T) */
   const double pressure = 1010.0; /* pressure (hPa) */
   const short int ref_option = 2; /*0 ... no refraction, 1 ... include refraction, using 'standard' atmospheric
    conditions, 2 ... include refraction, using atmospheric parameters input in the 'location' structure. */

    //observer on the near-Earth spacecraft1. Both input vectors are with respect to true equator and equinox of date.
   double rsun=696000.0, ddrsun=0.0;
   double rearth=6371.0084; //km
   double he=height/1000.0;
   double dtheta=0.0;
   double jd_beg, jd_end, jd_utc, jd_tt, jd_ut1, x, secdif, jd_tdb, delta_t, zd, az, rar, decr;
   
   //structures
   observer obs_geoc, obs_loc ;
   cat_entry cat_star, dummy_star;
   //object star, moon, mars, neptune, sun,venus;
   object  moon, mars, neptune, sun, venus;
   sky_pos t_place;


/* Make a structure of the observer */
   make_observer_at_geocenter ( &obs_geoc );
   make_observer_on_surface (latitude, longitude, height, temperature, pressure, &obs_loc);
    
/* Make structures of type 'object' for the Star, Sun, Moon, Mars and earth. */
  //  make_cat_entry ("647080","SC10MA",647080,0.001823085,25.88645705, 20.23, -7.14, 4.34, -31.0, &cat_star);
  // make_cat_entry("2111805", "SC10MA",211805,0.0,0.0,0.0,0.0,0.0,0.0, &cat_star);
   make_cat_entry("2111805", "xxx", 2111805,  6.39919, - 52.69566146 ,    19.92,     23.21 ,  10.56 ,  20.50, &cat_star);
   make_cat_entry("DUMMY", "xxx", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &dummy_star);
   /* if ((error = make_object (2,0,"star",&cat_star, &star)) != 0)
    {
        printf ("Error %d from make_object (cat_star)\n", error);
        return (error);
    }*/
	if ((error = make_object(0, 2, "Venus", &dummy_star, &venus)) != 0) // Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11
	{
		printf("Error %d from make_object (Venus)\n", error);
		return (error);
	}
  if ((error = make_object (0,11,"Moon",&dummy_star, &moon)) != 0) // Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11
   {
      printf ("Error %d from make_object (Moon)\n", error);
      return (error);
   }
   if ((error = make_object (0,4,"Mars",&dummy_star, &mars)) != 0)
   {
      printf ("Error %d from make_object (Mars)\n", error);
      return (error);
   }
    if ((error = make_object (0,8,"Neptune",&dummy_star, &neptune)) != 0)
    {
        printf ("Error %d from make_object (Neptune)\n", error);
        return (error);
    }
    if ((error = make_object (0,10,"Sun",&dummy_star, &sun)) != 0)
    {
        printf ("Error %d from make_object (Sun)\n", error);
        return (error);
    }

/* Open the JPL binary ephemeris file, here named "JPLEPH". */
   if ((error = ephem_open ("JPLEPH", &jd_beg,&jd_end,&de_num)) != 0)
   {
      if (error == 1)
         printf ("JPL ephemeris file not found.\n");
       else
         printf ("Error reading JPL ephemeris file header.\n");
      return (error);
   }
    else
   {
     // printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
        // de_num, jd_beg, jd_end);
      //printf ("\n");
   }
   
/* Angle between the refracted ture horizon and the astronomical horizon */
    dtheta=sqrt(2.0*rearth*he+he*he)/rearth*(1.0-1.0/7.0)*RAD2DEG; //deg
    //printf(" dtheta=%f \n",dtheta);
//求一下大致的日出时间，方便设定迭代初始值
	double delt,b1;
	short int year1, month1, day1;
	double hour1,jd_utc1,dt,ha;//ha:时角
	year1 = year; month1 = 1; day1 = 1; hour1 = 0.0;//一年的开始时间
	jd_utc = julian_date(year, month, day, hour);
	jd_utc1 = julian_date(year1, month1, day1, hour1);
	dt = jd_utc - jd_utc1;
	//printf("dt:%10.9f\n", dt);
	b1 = TWOPI*(dt - 1) / 365.2422;
	delt = 0.006918 - 0.399912*cos(b1) + 0.070257*sin(b1) - 0.006758*cos(2 * b1) + 0.000907*sin(2 * b1) - 0.002697*cos(3 * b1) + 0.00148*sin(3 * b1);//太阳赤纬
	//printf("太阳赤纬:%10.9f\n", delt*RAD2DEG);
	ha = acos(-tan(latitude*DEG2RAD)*tan(delt));//球面三角公式
	//printf("时角:%10.9f\n",ha*RAD2DEG/15.0);
	hour = 12.0 - ha*RAD2DEG / 15.0;//午时12点是视太阳时，减去日出日落时角
	//printf("hour:%10.9f\n", hour-(longitude-120.0)/15.0-16.0/60.0);
	//printf("方位角:%10.9f\n", acos(sin(-17.4*DEG2RAD)/cos(latitude*DEG2RAD))*RAD2DEG);
	//hour -= longitude / 15.0;//求格林尼治时间
	hour -= 25.0 / 60.0;//考虑均时差
	jd_utc = julian_date(year, month, day, hour)+1;//求第二天的日出时间
	jd_utc -= longitude / 15.0/24.0;
	/* Establish time arguments. */
	do{
		jd_utc+= 0.1/86400;
		jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		jd_ut1 = jd_utc + ut1_utc / 86400.0;
		delta_t = 32.184 + leap_secs - ut1_utc;
		secdif = 0.0;
		tdb2tt(jd_tt, &x, &secdif);
		jd_tdb = jd_tt + secdif / 86400.0;
		double rat = 0.0, dect = 0.0;
		if ((error = place(jd_tt, &sun, &obs_loc, delta_t, coord_sys, accuracy,
			&t_place)) != 0)
		{
			printf("Error %d from place.", error);
			return (error);
		}
		ddrsun = rsun / (t_place.dis*AU_KM) / DEG2RAD;
		equ2hor(jd_ut1, delta_t, accuracy, x_pole, y_pole, &obs_loc.on_surf, t_place.ra, t_place.dec, ref_option, &zd, &az, &rar, &decr);
		//printf("高度角:%10.9f\n", 90.0-zd);
		//printf("dtheta:%10.9f\n", fabs(zd - 90.0 - ddrsun - dtheta));
	} while (fabs(zd - 90.0 - ddrsun - dtheta) > ddzz);
    //UTC to Beijing Time
    cal_date (jd_utc+8.0/24.0, &year, &month, &day, &hour);
	a = (int)hour;
	b = (int)((hour - a) * 60);
	c = (hour - a - b / 60.0)*3600.0;
	printf("地理坐标经纬度为:(%6.4f,%6.4f)的地方，北京日期: %4d-%02d-%02d号的\n",longitude,latitude,year, month, day);
	printf("(1)日出方位角  ：%6.4f\n", az);
	printf("(2)日出北京时间：%02d:%02d:%4.2f\n", a, b, c);
    ephem_close();
   return (0);
}
