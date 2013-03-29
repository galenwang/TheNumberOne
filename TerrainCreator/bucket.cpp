#include "stdafx.h"
#include "bucket.h"
#include "TerrainProcessingUtils.h"
#include "TexCoordPreserver.h"
using std::string;
FILE *ft;

extern std::string cultures[];
extern int culture_cc;
extern Geodesy geoCenter;
extern Carton catCenter;
extern double gbs_radius;
extern std::vector<SimPoint> CultPoints;
extern std::vector<SimPoint> TreePoints;
extern void AddCorrectionSpaceMesh(osg::ref_ptr<osg::Group> gp);
extern osg::Matrix coordMatrix;
extern osg::Matrix rotateMatrix;
extern int FindIdxByColor(osg::Vec4 color);
extern std::vector<ColorName> colors;
extern char tempDirectory[_MAX_PATH];
extern std::vector<std::string> allDetailedTxt;
extern char tempDirectory[_MAX_PATH];
extern char debugFile[_MAX_PATH];
extern double GetElevation(double lon, double lat, double pre);
extern osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToTranglesFace(osg::ref_ptr<osg::Vec3Array> gvx);
extern osg::ref_ptr<osg::Geometry> ConvertFromLINELOOPToSkirtFace(osg::ref_ptr<osg::Vec3Array> gvx);
extern osg::ref_ptr<osg::Group> allLakeSkirt;
extern osg::ref_ptr<osg::Group> allRoadSkirt;
extern osg::ref_ptr<osg::Group> allDynamicGrass;
extern osg::ref_ptr<osg::Group> allIslands;
extern float GetWidthFromName(std::string matname);
extern std::vector<std::string> nameParts;
extern osg::Vec3 getSpaceVertex(osg::Vec3 v);
extern void getProjectionMapByVertex(ProjectionMap *pm, osg::Vec3 v);

void AnalysisBTG(char *btgName);
void WriteBTGHeader();
void CreateRandomPointsFromGroup(osg::ref_ptr<osg::Group> frst, int mode);
float  mWoodCoverage;
std::vector<std::string> materialParts;
int splitMaterialName(string matname);
void InitMMLLEx(MaxMinLonLat *m, string fn);
void SetMMLLEx(MaxMinLonLat *m, double lon, double lat);

/******************************************************************
 *
 *    ��������ɻ����㷨
 *
 ******************************************************************/

mt random_seed;

#define MT(i) mt->array[i]

void mt_init(mt *mt, unsigned int seed)
{
    int i;
    MT(0)= seed;
    for(i=1; i<MT_N; i++)
        MT(i) = (1812433253 * (MT(i-1) ^ (MT(i-1) >> 30)) + i);
    mt->index = MT_N+1;
}

unsigned int mt_rand32(mt *mt)
{
    unsigned int i, y;
    if(mt->index >= MT_N) {
        for(i=0; i<MT_N; i++) {
            y = (MT(i) & 0x80000000) | (MT((i+1)%MT_N) & 0x7fffffff);
            MT(i) = MT((i+MT_M)%MT_N) ^ (y>>1) ^ (y&1 ? 0x9908b0df : 0);
        }
        mt->index = 0;
    }
    y = MT(mt->index++);
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);
    return y;
}

double mt_rand(mt *mt)
{
    /* divided by 2^32-1 */ 
    return (double)mt_rand32(mt) * (1.0/4294967295.0); 
}



/******************************************************************
 *
 *    ����ϵת�������롢�ǶȵȵȻ����㷨����
 *
 ******************************************************************/

/* �Ƕ�ת���� */
double deg2rad(double val)
{
	return val*SGD_DEGREES_TO_RADIANS;
}

/* ����ת�Ƕ� */
double rad2deg(double val)
{
	return val*SGD_RADIANS_TO_DEGREES;
}

/* ����Բ�����ɵ�ֵ */
double pi()
{
	return SGD_PI;
}

/* ����lat1, lon1, lat2, lon2, ���㿪ʼ��ͽ�����;��롣 */
// given lat1, lon1, lat2, lon2, calculate starting and ending
// az1, az2 and distance (s).  Lat, lon, and azimuth are in degrees.
// distance in meters
static int _geo_inverse_wgs_84( double lat1, double lon1, double lat2,
			double lon2, double *az1, double *az2,
                        double *s )
{
    double a = _EQURAD, rf = _FLATTENING;
    int iter=0;
    double testv = 1.0E-10;
    double f = ( rf > 0.0 ? 1.0/rf : 0.0 );
    double b = a*(1.0-f);
    // double e2 = f*(2.0-f); // unused in this routine
    double phi1 = deg2rad(lat1), lam1 = deg2rad(lon1);
    double sinphi1 = sin(phi1), cosphi1 = cos(phi1);
    double phi2 = deg2rad(lat2), lam2 = deg2rad(lon2);
    double sinphi2 = sin(phi2), cosphi2 = cos(phi2);
	
    if( (fabs(lat1-lat2) < testv && 
	 ( fabs(lon1-lon2) < testv)) || (fabs(lat1-90.0) < testv ) )
    {	
	// TWO STATIONS ARE IDENTICAL : SET DISTANCE & AZIMUTHS TO ZERO */
	*az1 = 0.0; *az2 = 0.0; *s = 0.0;
	return 0;
    } else if(  fabs(cosphi1) < testv ) {
	// initial point is polar
	int k = _geo_inverse_wgs_84( lat2,lon2,lat1,lon1, az1,az2,s );
	k = k; // avoid compiler error since return result is unused
	b = *az1; *az1 = *az2; *az2 = b;
	return 0;
    } else if( fabs(cosphi2) < testv ) {
	// terminal point is polar
        double _lon1 = lon1 + 180.0f;
	int k = _geo_inverse_wgs_84( lat1, lon1, lat1, _lon1, 
				    az1, az2, s );
	k = k; // avoid compiler error since return result is unused
	*s /= 2.0;
	*az2 = *az1 + 180.0;
	if( *az2 > 360.0 ) *az2 -= 360.0; 
	return 0;
    } else if( (fabs( fabs(lon1-lon2) - 180 ) < testv) && 
	       (fabs(lat1+lat2) < testv) ) 
    {
	// Geodesic passes through the pole (antipodal)
	double s1,s2;
	_geo_inverse_wgs_84( lat1,lon1, lat1,lon2, az1,az2, &s1 );
	_geo_inverse_wgs_84( lat2,lon2, lat1,lon2, az1,az2, &s2 );
	*az2 = *az1;
	*s = s1 + s2;
	return 0;
    } else {
	// antipodal and polar points don't get here
	double dlam = lam2 - lam1, dlams = dlam;
	double sdlams,cdlams, sig,sinsig,cossig, sinaz,
	    cos2saz, c2sigm;
	double tc,temp, us,rnumer,denom, ta,tb;
	double cosu1,sinu1, sinu2,cosu2;

	// Reduced latitudes
	temp = (1.0-f)*sinphi1/cosphi1;
	cosu1 = 1.0/sqrt(1.0+temp*temp);
	sinu1 = temp*cosu1;
	temp = (1.0-f)*sinphi2/cosphi2;
	cosu2 = 1.0/sqrt(1.0+temp*temp);
	sinu2 = temp*cosu2;
    
	do {
	    sdlams = sin(dlams), cdlams = cos(dlams);
	    sinsig = sqrt(cosu2*cosu2*sdlams*sdlams+
			  (cosu1*sinu2-sinu1*cosu2*cdlams)*
			  (cosu1*sinu2-sinu1*cosu2*cdlams));
	    cossig = sinu1*sinu2+cosu1*cosu2*cdlams;
	    
	    sig = atan2(sinsig,cossig);
	    sinaz = cosu1*cosu2*sdlams/sinsig;
	    cos2saz = 1.0-sinaz*sinaz;
	    c2sigm = (sinu1 == 0.0 || sinu2 == 0.0 ? cossig : 
		      cossig-2.0*sinu1*sinu2/cos2saz);
	    tc = f*cos2saz*(4.0+f*(4.0-3.0*cos2saz))/16.0;
	    temp = dlams;
	    dlams = dlam+(1.0-tc)*f*sinaz*
		(sig+tc*sinsig*
		 (c2sigm+tc*cossig*(-1.0+2.0*c2sigm*c2sigm)));
	    if (fabs(dlams) > pi() && iter++ > 50) {
		return iter;
	    }
	} while ( fabs(temp-dlams) > testv);

	us = cos2saz*(a*a-b*b)/(b*b); // !!
	// BACK AZIMUTH FROM NORTH
	rnumer = -(cosu1*sdlams);
	denom = sinu1*cosu2-cosu1*sinu2*cdlams;
	*az2 = rad2deg(atan2(rnumer,denom));
	if( fabs(*az2) < testv ) *az2 = 0.0;
	if(*az2 < 0.0) *az2 += 360.0;

	// FORWARD AZIMUTH FROM NORTH
	rnumer = cosu2*sdlams;
	denom = cosu1*sinu2-sinu1*cosu2*cdlams;
	*az1 = rad2deg(atan2(rnumer,denom));
	if( fabs(*az1) < testv ) *az1 = 0.0;
	if(*az1 < 0.0) *az1 += 360.0;

	// Terms a & b
	ta = 1.0+us*(4096.0+us*(-768.0+us*(320.0-175.0*us)))/
	    16384.0;
	tb = us*(256.0+us*(-128.0+us*(74.0-47.0*us)))/1024.0;

	// GEODETIC DISTANCE
	*s = b*ta*(sig-tb*sinsig*
		   (c2sigm+tb*(cossig*(-1.0+2.0*c2sigm*c2sigm)-tb*
			       c2sigm*(-3.0+4.0*sinsig*sinsig)*
			       (-3.0+4.0*c2sigm*c2sigm)/6.0)/
		    4.0));
	return 0;
    }
}

/* ����������γ������ϵ��֮��ľ��룬��λΪ�� */
double distanceM(Geodesy *p1, Geodesy *p2)
{
  double course1, course2, distance;
  int r = _geo_inverse_wgs_84(p1->_lat, p1->_lon,
                                p2->_lat, p2->_lon,
                                &course1, &course2, &distance);
  if (r != 0) {
    return 0.0;
  }
  
  return distance;
}

/* ����ֱ������ϵת��Ϊ��γ������ϵ */
void CartToGeod(Carton *cart, Geodesy *geod)
{
  double X, Y, Z, XXpYY, sqrtXXpYY, p, q, r, s, t, u, v, w, k, D, sqrtDDpZZ;	
  double a = _EQURAD;
  double ra2 = 1/(_EQURAD*_EQURAD);
  double e2 = E2;
  double e4 = E2*E2;
	
  // according to
  // H. Vermeille,
  // Direct transformation from geocentric to geodetic ccordinates,
  // Journal of Geodesy (2002) 76:451-454
  X = cart->x;
  Y = cart->y;
  Z = cart->z;
  XXpYY = X*X+Y*Y;
  if( XXpYY + Z*Z < 25 ) {
    // This function fails near the geocenter region, so catch that special case here.
    // Define the innermost sphere of small radius as earth center and return the 
    // coordinates 0/0/-EQURAD. It may be any other place on geoide's surface,
    // the Northpole, Hawaii or Wentorf. This one was easy to code ;-)
    geod->_lon = 0.0 ;
    geod->_lat = 0.0 ;
    geod->_elevation = -_EQURAD ;
  }
    
  sqrtXXpYY = sqrt(XXpYY);
  p = XXpYY*ra2;
  q = Z*Z*(1-e2)*ra2;
  r = 1/6.0*(p+q-e4);
  s = e4*p*q/(4*r*r*r);
/* 
  s*(2+s) is negative for s = [-2..0]
  slightly negative values for s due to floating point rounding errors
  cause nan for sqrt(s*(2+s))
  We can probably clamp the resulting parable to positive numbers
*/
  if( s >= -2.0 && s <= 0.0 )
    s = 0.0;
  t = pow(1+s+sqrt(s*(2+s)), 1/3.0);
  u = r*(1+t+1/t);
  v = sqrt(u*u+e4*q);
  w = e2*(u+v-q)/(2*v);
  k = sqrt(u+v+w*w)-w;
  D = k*sqrtXXpYY/(k+e2);
  geod->_lon = (2*atan2(Y, X+sqrtXXpYY))*SGD_RADIANS_TO_DEGREES;
  sqrtDDpZZ = sqrt(D*D+Z*Z);
  geod->_lat = (2*atan2(Z, D+sqrtDDpZZ))*SGD_RADIANS_TO_DEGREES;
  geod->_elevation = ((k+e2-1)*sqrtDDpZZ/k);
}

/* ��γ������ϵת��Ϊ����ֱ������ϵ */
void GeodToCart(Geodesy *geod, Carton *cart)
{
  double a = _EQURAD;
  double e2 = E2;
  // according to
  // H. Vermeille,
  // Direct transformation from geocentric to geodetic ccordinates,
  // Journal of Geodesy (2002) 76:451-454
  double lambda = deg2rad(geod->_lon);
  double phi = deg2rad(geod->_lat);
  double h = geod->_elevation;
  double sphi = sin(phi);
  double n = a/sqrt(1-e2*sphi*sphi);
  double cphi = cos(phi);
  double slambda = sin(lambda);
  double clambda = cos(lambda);
  cart->x = (h+n)*cphi*clambda;
  cart->y = (h+n)*cphi*slambda;
  cart->z = (h+n-e2*n)*sphi;
}

/* ��������ֱ������ϵ��֮��ľ���ƽ���� */
double distance3Dsquared(Carton *a, Carton *b)
{
    double x, y, z;

    x = a->x - b->x;
    y = a->y - b->y;
    z = a->z - b->z;

    return(x*x + y*y + z*z);
}

/* ��������ֱ������ϵ��֮��ľ��� */
double distance3D(Carton *a, Carton *b)
{
    double x, y, z;

    x = a->x - b->x;
    y = a->y - b->y;
    z = a->z - b->z;

    return sqrt(x*x + y*y + z*z);
}


/******************************************************************
 *
 *    ���ι���Ԫ Bucket �Ļ����㷨��Ϊ�˱��ڶԵؾ�ϸ�ڵĹ�����һ����γ�Ȼ��ֳ����ɸ�bucket��
 *
 ******************************************************************/

// return the horizontal tile span factor based on latitude
static double sg_bucket_span( double l ) {
    if ( l >= 89.0 ) {
	return 360.0;
    } else if ( l >= 88.0 ) {
	return 8.0;
    } else if ( l >= 86.0 ) {
	return 4.0;
    } else if ( l >= 83.0 ) {
	return 2.0;
    } else if ( l >= 76.0 ) {
	return 1.0;
    } else if ( l >= 62.0 ) {
	return 0.5;
    } else if ( l >= 22.0 ) {
	return 0.25;
    } else if ( l >= -22.0 ) {
	return 0.125;
    } else if ( l >= -62.0 ) {
	return 0.25;
    } else if ( l >= -76.0 ) {
	return 0.5;
    } else if ( l >= -83.0 ) {
	return 1.0;
    } else if ( l >= -86.0 ) {
	return 2.0;
    } else if ( l >= -88.0 ) {
	return 4.0;
    } else if ( l >= -89.0 ) {
	return 8.0;
    } else {
	return 360.0;
    }
}

/* ���ɵ�ǰbucket����ֵ */
long int gen_index(Bucket *B) 
{
	return ((B->lon + 180) << 14) + ((B->lat + 90) << 6) + (B->y << 3) + B->x;
}

/* ���ɵ�ǰbucket����ֵ�ַ��� */
void gen_index_str(Bucket *B, char *idx) 
{
	sprintf(idx, "%d", ((B->lon + 180) << 14) + ((B->lat + 90) << 6) + (B->y << 3) + B->x);
}

/* ��ȡ��ǰbucket���ĵ�γ��ֵ */
double get_center_lat(Bucket *B) 
{
    return B->lat + B->y / 8.0 + SG_HALF_BUCKET_SPAN;
}

// return width of the tile in degrees
double get_width(Bucket *B) 
{
    return sg_bucket_span( get_center_lat(B) );
}


// return height of the tile in degrees
double get_height(Bucket *B) 
{
    return SG_BUCKET_SPAN;
}

// Build the path name for this bucket
void gen_base_path(Bucket *B, char *path)
{
    // long int index;
    int top_lon, top_lat, main_lon, main_lat;
    char hem, pole;
    char raw_path[256];

    top_lon = B->lon / 10;
    main_lon = B->lon;
    if ( (B->lon < 0) && (top_lon * 10 != B->lon) ) {
		top_lon -= 1;
    }
    top_lon *= 10;
    if ( top_lon >= 0 ) {
		hem = 'e';
    } else {
		hem = 'w';
		top_lon *= -1;
    }
    if ( main_lon < 0 ) {
		main_lon *= -1;
    }
    
    top_lat = B->lat / 10;
    main_lat = B->lat;
    if ( (B->lat < 0) && (top_lat * 10 != B->lat) ) {
		top_lat -= 1;
    }
    top_lat *= 10;
    if ( top_lat >= 0 ) {
		pole = 'n';
    } else {
		pole = 's';
		top_lat *= -1;
    }
    if ( main_lat < 0 ) {
		main_lat *= -1;
    }

    sprintf(raw_path, "%c%03d%c%02d\\%c%03d%c%02d", 
	    hem, top_lon, pole, top_lat, 
	    hem, main_lon, pole, main_lat);

    strcpy(path, raw_path);
}

/* ��ȡ��ǰbucket���ĵ㾭��ֵ */
double get_center_lon(Bucket *B)
{
   double span = sg_bucket_span( B->lat + B->y / 8.0 + SG_HALF_BUCKET_SPAN );

   if ( span >= 1.0 ) {
      return B->lon + span / 2.0;
   } else {
      return B->lon + B->x * span + span / 2.0;
   }
}

/* ��������һ��ľ�γ��ֵ����bucket */
void set_bucket(Bucket *B, double dlon, double dlat ) 
{
    //
    // latitude first
    //
    double span = sg_bucket_span( dlat );
    double diff = dlon - (double)(int)dlon;

    // cout << "diff = " << diff << "  span = " << span << endl;
    if ( (dlon >= 0) || (fabs(diff) < SG_EPSILON) ) {
	B->lon = (int)dlon;
    } else {
	B->lon = (int)dlon - 1;
    }

    // find subdivision or super lon if needed
    if ( span < SG_EPSILON ) {
	// polar cap
	B->lon = 0;
	B->x = 0;
    } else if ( span <= 1.0 ) {
	B->x = (int)((dlon - B->lon) / span);
    } else {
	if ( dlon >= 0 ) {
	    B->lon = (int)( (int)(B->lon / span) * span);
	} else {
	    // cout << " lon = " << lon 
	    //  << "  tmp = " << (int)((lon-1) / span) << endl;
	    B->lon = (int)( (int)((B->lon + 1) / span) * span - span);
	    if ( B->lon < -180 ) {
		B->lon = -180;
	    }
	}
	B->x = 0;
    }

    //
    // then latitude
    //
    diff = dlat - (double)(int)dlat;

    if ( (dlat >= 0) || (fabs(diff) < SG_EPSILON) ) {
	B->lat = (int)dlat;
    } else {
	B->lat = (int)dlat - 1;
    }
    B->y = (int)((dlat - B->lat) * 8);
}

/* ��������һ������ֵ����bucket */
void SGBucket(Bucket *B, const long int bindex) 
{
    long int index = bindex;
	
    B->lon = index >> 14;
    index -= B->lon << 14;
    B->lon -= 180;

    B->lat = index >> 6;
    index -= B->lat << 6;
    B->lat -= 90;

    B->y = index >> 3;
    index -= B->y << 3;

    B->x = index;
}


/**********************************************************
 *
 *  
 *      BTG ��һЩ���̣�btg�ǵؾ��ļ����ͣ�һ��ר�����͵���ά�����ļ����ɵ����湹�ɡ�
 *
 *
 **********************************************************/

std::vector<TerPoint> curPnt;					//�����������洢��ά�ؾ����еĵ㡣
std::vector<TerTextcoord> curText;				//�����������洢��ά�������е���������ꡣ
std::vector<TriType> curTri;					//�����������洢��ά�������е��档
std::vector<string> myCultures;	

/* btg�ļ���һЩ��ر��� */
double dXpart, dYpart, dZpart;
float  fRadiusofBS;
int iCreationTime, totalNormal, numNew, totalObjects;
unsigned short iVersion, iMagicNum, iNumOfToplvlObject; 
double gbCenterX, gbCenterY, gbCenterZ;

int WriteDataStructToBTG_OnCommand(char *fname);

/* ��ȡbtg�ļ����洢��������ݽṹ�����С� */
int ReadBTGToDataStruct(char *fname)
{
	int m, i, j, k, indexType;
	char meterial[32];
	gzFile fp;

	//��ʼ���ڴ�
	curPnt.clear();
	curText.clear();
	curTri.clear();
	totalNormal = numNew = totalObjects = 0;
	memset(meterial, 0, sizeof(meterial));

	//���ļ�����ȡ
	fp = gzopen(fname, "rb");
	if(fp == NULL) return 0;

	//Head
	gzread(fp, &iVersion, sizeof(short));
	gzread(fp, &iMagicNum, sizeof(short));
	gzread(fp, &iCreationTime, sizeof(int));
	gzread(fp, &iNumOfToplvlObject, sizeof(short));

	for(m=0;m<iNumOfToplvlObject;m++)
	{
		//Object Head
		char cObjType;
		unsigned short iNumOfObjProps, iNumOfObjElems;
		int iLength;

		gzread(fp, &cObjType, sizeof(char));
		gzread(fp, &iNumOfObjProps, sizeof(short));
		gzread(fp, &iNumOfObjElems, sizeof(short));

		if(iNumOfObjProps)
		{
			char cPropType;
			char sProps[64];
			int iPropLength;

			memset(sProps, 0, 64);
			for(i=0;i<iNumOfObjProps;i++)
			{
				gzread(fp, &cPropType, sizeof(char));
				gzread(fp, &iPropLength, sizeof(int));
				switch(cPropType)
				{
				case 0:		//Material
					{
						gzread(fp, sProps, sizeof(char)*iPropLength);
						strcpy(meterial, sProps);
					}
					break;
				case 1:		//Index Types
					{
						gzread(fp, sProps, sizeof(char)*iPropLength);
						indexType = sProps[0];
					}
					break;
				}
			}
		}

		TriType tt;
		tt.mElement.clear();
		for(k=0;k<iNumOfObjElems;k++)
		{
			gzread(fp, &iLength, sizeof(int));

			switch(cObjType)
			{
			case 0:		// Bounding Sphere  
				{
					gzread(fp, &dXpart, sizeof(double));
					gzread(fp, &dYpart, sizeof(double));
					gzread(fp, &dZpart, sizeof(double));
					gzread(fp, &fRadiusofBS, sizeof(float));

					gbCenterX = dXpart; gbCenterY = dYpart; gbCenterZ = dZpart;
				}
				break;
			case 1:		// Vertex List  
				{
					int iNumberOfVertex = iLength / 12;
					float fVx, fVy, fVz;
					Carton cart;
					Geodesy geod;

					j=0;
					for(i=0;i<iNumberOfVertex;i++)
					{
						gzread(fp, &fVx, sizeof(float));
						gzread(fp, &fVy, sizeof(float));
						gzread(fp, &fVz, sizeof(float));

						cart.x = fVx + gbCenterX; cart.y = fVy + gbCenterY; cart.z = fVz + gbCenterZ;
						CartToGeod(&cart, &geod);

						TerPoint tp;
						tp.px   = cart.x; tp.py   = cart.y; tp.pz   = cart.z; tp.lon  = geod._lon; tp.lat  = geod._lat; tp.elev = geod._elevation; 
						curPnt.push_back(tp);
					}
				}	
				break;
			case 2:		// Normal List  
				{
					int iNumberOfNormal = iLength / 3;
					char bNx, bNy, bNz;

					j=0;
					for(i=0;i<iNumberOfNormal;i++)
					{
						gzread(fp, &bNx, sizeof(char));
						gzread(fp, &bNy, sizeof(char));
						gzread(fp, &bNz, sizeof(char));

						curPnt[totalNormal].nx = bNx; 	curPnt[totalNormal].ny = bNy; curPnt[totalNormal].nz = bNz; 
						totalNormal++;
					}
				}	
				break;
			case 3:		// Texture Coordinates List  
				{
					int iNumberOfTexture = iLength / 8;
					float fTx, fTy;

					j=0;
					for(i=0;i<iNumberOfTexture;i++)
					{
						gzread(fp, &fTx, sizeof(float));
						gzread(fp, &fTy, sizeof(float));

						TerTextcoord tt; tt.fU = fTx; tt.fV = fTy; curText.push_back(tt);
					}
				}	
				break;
			case 4:		// Color List  
				{
					int iNumberOfColor = iLength / 16;
					float fVR, fVG, fVB, fVA;

					j=0;
					for(i=0;i<iNumberOfColor;i++)
					{
						gzread(fp, &fVR, sizeof(float));
						gzread(fp, &fVG, sizeof(float));
						gzread(fp, &fVB, sizeof(float));
						gzread(fp, &fVA, sizeof(float));
					}
				}	
				break;
			case 9:		// Points  d
			case 10:	// Individual Triangles  
			case 11:	// Triangle Strips  
			case 12:	// Triangle Fans  
				{
					int iNumberOfPoints;
					unsigned short sPointV, sPointN, sPointT;

					TriElement te; te.mVertex.clear(); te.minP = 32767;
					if((indexType == 3)||(indexType == 5)||(indexType == 9))
					{
						iNumberOfPoints = iLength / 4;

						j=0;
						for(i=0;i<iNumberOfPoints;i++)
						{
							gzread(fp, &sPointV, sizeof(short));
							gzread(fp, &sPointT, sizeof(short));

							TriVertex tv;
							tv.P = sPointV; tv.N = 0; tv.T = sPointT; te.mVertex.push_back(tv);
							if(sPointV < te.minP) te.minP = sPointV;
						}
					}

					if((indexType == 7)||(indexType == 11)||(indexType == 13))
					{
						iNumberOfPoints = iLength / 6;

						j=0;
						for(i=0;i<iNumberOfPoints;i++)
						{
							gzread(fp, &sPointV, sizeof(short));
							gzread(fp, &sPointN, sizeof(short));
							gzread(fp, &sPointT, sizeof(short));

							TriVertex tv; tv.P = sPointV; tv.N = sPointN; tv.T = sPointT; te.mVertex.push_back(tv);
							if(sPointV < te.minP) te.minP = sPointV;
						}
					}
					te.iSum = te.mVertex.size();
					tt.mElement.push_back(te);
				}	
				break;
			}// switch
		}// for k
		
		if((cObjType == 9)||(cObjType == 10)||(cObjType == 11)||(cObjType == 12))
		{
			totalObjects += tt.mElement.size();
			tt.iElems = tt.mElement.size();
			strcpy(tt.cMtype, meterial);
			tt.cObjType = cObjType;
			tt.cPropVal = indexType;
			curTri.push_back(tt);
		}		
	}// for m 

	gzclose(fp);
	return 1;
}

/* ��������ݽṹ�����е�����д�뵽btg�ļ����� */
int WriteDataStructToBTG(char *fname)
{
	int m, k, indexType;
	char cPropType, cPropVal;
	int iPropLength;
	char cObjType;
	unsigned short iNumOfObjProps, iNumOfObjElems;
	int iLength;
	gzFile fo;

	fo = gzopen(fname, "wb");
	if(fo == NULL) return 0;

	/* дbtg�ļ�ͷ */
	gzwrite(fo, &iVersion, sizeof(short));
	gzwrite(fo, &iMagicNum, sizeof(short));
	gzwrite(fo, &iCreationTime, sizeof(int));
	gzwrite(fo, &iNumOfToplvlObject, sizeof(short));

	/* дbtg�ļ���Χ�� */
	cObjType = 0; iNumOfObjProps = 0; iNumOfObjElems = 1; iLength = 0x1C;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{

		dXpart = gbCenterX;
		dYpart = gbCenterY;
		dZpart = gbCenterZ;

		gzwrite(fo, &dXpart, sizeof(double));
		gzwrite(fo, &dYpart, sizeof(double));
		gzwrite(fo, &dZpart, sizeof(double));
		gzwrite(fo, &fRadiusofBS, sizeof(float));
	}

	/* дbtg�ļ����ж��㣬��¼��Χ�����ĵ���������Ĳ�ֵ */
	cObjType = 1; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curPnt.size() * 12;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fVx, fVy, fVz;

		for(m=0;m<(int)curPnt.size();m++)
		{
			fVx = curPnt[m].px - gbCenterX;
			fVy = curPnt[m].py - gbCenterY;
			fVz = curPnt[m].pz - gbCenterZ;

			gzwrite(fo, &fVx, sizeof(float));
			gzwrite(fo, &fVy, sizeof(float));
			gzwrite(fo, &fVz, sizeof(float));
		}
	}

	/* дbtg�ļ������������ɫֵ */
	cObjType = 4; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = 0;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	
	/* дbtg�ļ���������ķ�����ֵ */
	cObjType = 2; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = totalNormal * 3;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		char bNx, bNy, bNz;

		for(m=0;m<totalNormal;m++)
		{
			bNx = curPnt[m].nx; 
			bNy = curPnt[m].ny; 
			bNz = curPnt[m].nz;
			
			gzwrite(fo, &bNx, sizeof(char));
			gzwrite(fo, &bNy, sizeof(char));
			gzwrite(fo, &bNz, sizeof(char));
		}
	}

	/* дbtg�ļ������������������ֵ */
	cObjType = 3; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curText.size() * 8;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fTx, fTy;
		
		for(m=0;m<(int)curText.size();m++)
		{
			fTx = curText[m].fU; 
			fTy = curText[m].fV;
			
			gzwrite(fo, &fTx, sizeof(float));
			gzwrite(fo, &fTy, sizeof(float));
		}
	}

	/* дbtg�ļ��������棬������������ɣ�type:10Ϊ���������Σ�11Ϊ���Σ�12Ϊ���� */
	for(unsigned int i=0;i<curTri.size();i++)
	{
		cObjType = curTri[i].cObjType; iNumOfObjProps = 2; iNumOfObjElems = curTri[i].iElems; 
		gzwrite(fo, &cObjType, sizeof(char));
		gzwrite(fo, &iNumOfObjProps, sizeof(short));
		gzwrite(fo, &iNumOfObjElems, sizeof(short));

		//for(j=0;j<strlen(curTri[i].cMtype);j++) { if(curTri[i].cMtype[j]>0x60) curTri[i].cMtype[j]-=0x20; }		//�����������е�Сд��ĸ�ĳɴ�д��
		cPropType = 0; iPropLength = strlen(curTri[i].cMtype);
		gzwrite(fo, &cPropType, sizeof(char));
		gzwrite(fo, &iPropLength, sizeof(int));
		for(k=0;k<iPropLength;k++)  gzwrite(fo, &curTri[i].cMtype[k], sizeof(char));
		cPropType = 1; iPropLength = 1;  indexType = cPropVal = curTri[i].cPropVal;
		gzwrite(fo, &cPropType, sizeof(char));
		gzwrite(fo, &iPropLength, sizeof(int));
		gzwrite(fo, &cPropVal, sizeof(char));

		for(unsigned int j=0;j<curTri[i].mElement.size();j++)
		{
			unsigned short sPointV, sPointN, sPointT;

			if((indexType == 3)||(indexType == 5)||(indexType == 9))
			{
				iLength = curTri[i].mElement[j].mVertex.size() * 4;
				gzwrite(fo, &iLength, sizeof(int));

				for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size();k++)
				{
					sPointV = curTri[i].mElement[j].mVertex[k].P;
					if(indexType == 3) sPointT = curTri[i].mElement[j].mVertex[k].N;
					if(indexType == 9) sPointT = curTri[i].mElement[j].mVertex[k].T;

					gzwrite(fo, &sPointV, sizeof(short));
					gzwrite(fo, &sPointT, sizeof(short));
				}//k
			}

			if((indexType == 7)||(indexType == 11)||(indexType == 13))
			{
				iLength = curTri[i].mElement[j].mVertex.size() * 6;
				gzwrite(fo, &iLength, sizeof(int));

				for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size();k++)
				{
					sPointV = curTri[i].mElement[j].mVertex[k].P;
					sPointN = curTri[i].mElement[j].mVertex[k].N;
					sPointT = curTri[i].mElement[j].mVertex[k].T;

					gzwrite(fo, &sPointV, sizeof(short));
					gzwrite(fo, &sPointN, sizeof(short));
					gzwrite(fo, &sPointT, sizeof(short));
				}//k
			}
		}// j
	}// i

	gzclose(fo);
	return 1;
}

/* д��ͬһ������ˮ */
int WriteDataStructToBTG_River(char *fname)
{
	int m, k, indexType;
	char cPropType, cPropVal;
	int iPropLength;
	char cObjType;
	unsigned short iNumOfObjProps, iNumOfObjElems;
	int iLength;
	gzFile fo;

	fo = gzopen(fname, "wb");
	if(fo == NULL) return 0;

	iNumOfToplvlObject = 4;
	gzwrite(fo, &iVersion, sizeof(short));
	gzwrite(fo, &iMagicNum, sizeof(short));
	gzwrite(fo, &iCreationTime, sizeof(int));
	gzwrite(fo, &iNumOfToplvlObject, sizeof(short));

	cObjType = 0; iNumOfObjProps = 0; iNumOfObjElems = 1; iLength = 0x1C;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{

		dXpart = gbCenterX;
		dYpart = gbCenterY;
		dZpart = gbCenterZ;

		gzwrite(fo, &dXpart, sizeof(double));
		gzwrite(fo, &dYpart, sizeof(double));
		gzwrite(fo, &dZpart, sizeof(double));
		gzwrite(fo, &fRadiusofBS, sizeof(float));
	}

	cObjType = 1; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curPnt.size() * 12;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fVx, fVy, fVz;

		for(m=0;m<(int)curPnt.size();m++)
		{
			fVx = curPnt[m].px - gbCenterX;
			fVy = curPnt[m].py - gbCenterY;
			fVz = curPnt[m].pz - gbCenterZ;

			gzwrite(fo, &fVx, sizeof(float));
			gzwrite(fo, &fVy, sizeof(float));
			gzwrite(fo, &fVz, sizeof(float));
		}
	}

	cObjType = 3; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curText.size() * 8;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fTx, fTy;
		
		for(m=0;m<(int)curText.size();m++)
		{
			fTx = curText[m].fU; 
			fTy = curText[m].fV;
			
			gzwrite(fo, &fTx, sizeof(float));
			gzwrite(fo, &fTy, sizeof(float));
		}
	}

	iNumOfObjElems = 0;
	for(unsigned int i=0;i<curTri.size();i++)
	{
		if(curTri[i].cObjType != 9) iNumOfObjElems += curTri[i].iElems;
	}

	cObjType = 10; iNumOfObjProps = 2; 
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));

	char cMtype[10];
	sprintf(cMtype, "Stream");
	cPropType = 0; iPropLength = strlen(cMtype);
	gzwrite(fo, &cPropType, sizeof(char));
	gzwrite(fo, &iPropLength, sizeof(int));
	for(k=0;k<iPropLength;k++)  gzwrite(fo, &cMtype[k], sizeof(char));
	cPropType = 1; iPropLength = 1;  indexType = cPropVal = 9;
	gzwrite(fo, &cPropType, sizeof(char));
	gzwrite(fo, &iPropLength, sizeof(int));
	gzwrite(fo, &cPropVal, sizeof(char));
	
	for(unsigned int i=0;i<curTri.size();i++)
	{
		if(curTri[i].cObjType == 9) continue;
		for(unsigned int j=0;j<curTri[i].mElement.size();j++)
		{
			unsigned short sPointV, sPointT;

			iLength = curTri[i].mElement[j].mVertex.size() * 4;
			gzwrite(fo, &iLength, sizeof(int));

			for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size();k++)
			{
				sPointV = curTri[i].mElement[j].mVertex[k].P;
				sPointT = curTri[i].mElement[j].mVertex[k].T;

				gzwrite(fo, &sPointV, sizeof(short));
				gzwrite(fo, &sPointT, sizeof(short));
			}//k
		}// j
	}// i

	gzclose(fo);
	return 1;
}

/* ֻ���Ļ���Ϣ���ǲ���д��һ���µ�Btg�ļ� */
int WriteDataStructToBTG_OnlyCulture(char *fname)
{
	int m;
	/* ��cultures�����һ����Island�����Ա�Ѻ��ĵ����ڳ��� */
	myCultures.clear();
	for(m = 0; m < culture_cc; m++)
	{
		if(cultures[m] != "DynamicGrass")
			myCultures.push_back(cultures[m]);
	}

	/* ��ʼд�ļ� */
	if(!WriteDataStructToBTG_OnCommand(fname)) return 0;
	return 1;
}

/* ���Ļ���Ϣ��ϸ��������ǲ���д��һ���µ�Btg�ļ� */
int WriteDataStructToBTG_CulturePlusIsland(char *fname)
{
	int m;
	/* ��cultures�����ϸ�������һ����Island�����Ա�Ѻ��ĵ����ڳ��� */
	myCultures.clear();
	for(m = 0; m < culture_cc; m++)
		myCultures.push_back(cultures[m]);
	if(allDetailedTxt.size() > 0)
	{
		for(m = 0; m < (int)allDetailedTxt.size(); m++)
			myCultures.push_back(allDetailedTxt[m]);
	}
	myCultures.push_back("Island");	

	/* ��ʼд�ļ� */
	if(!WriteDataStructToBTG_OnCommand(fname)) return 0;
	return 1;
}

int WriteDataStructToBTG_OnCommand(char *fname)
{
	int m, j, k, indexType, cultNum = 0;
	char cPropType, cPropVal;
	int iPropLength;
	char cObjType;
	unsigned short iNumOfObjProps, iNumOfObjElems;
	int iLength;
	gzFile fo;
	
	/* ͳ���Ļ���ϢKind�ĸ��� */
	for(unsigned int i=0;i<curTri.size();i++)
	{
		for(j=0; j<(int)myCultures.size(); j++)
		{
			std::string curMtype = curTri[i].cMtype;
			if((!strcmp(curTri[i].cMtype, myCultures[j].c_str())) || (curMtype.substr(0, 7) == DEFINE_AIRPORT_STRING))
			{
				cultNum++;	break;
			}
		}
	}
	if(cultNum == 0) return 0;

	/* ��ʼдBTG�ļ� */
	fo = gzopen(fname, "wb");
	if(fo == NULL) return 0;
	
	iNumOfToplvlObject = 5 + cultNum;
	
	gzwrite(fo, &iVersion, sizeof(short));
	gzwrite(fo, &iMagicNum, sizeof(short));
	gzwrite(fo, &iCreationTime, sizeof(int));
	gzwrite(fo, &iNumOfToplvlObject, sizeof(short));

	cObjType = 0; iNumOfObjProps = 0; iNumOfObjElems = 1; iLength = 0x1C;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{

		dXpart = gbCenterX;
		dYpart = gbCenterY;
		dZpart = gbCenterZ;

		gzwrite(fo, &dXpart, sizeof(double));
		gzwrite(fo, &dYpart, sizeof(double));
		gzwrite(fo, &dZpart, sizeof(double));
		gzwrite(fo, &fRadiusofBS, sizeof(float));
	}

	cObjType = 1; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curPnt.size() * 12;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fVx, fVy, fVz;

		for(m=0;m<(int)curPnt.size();m++)
		{
			fVx = curPnt[m].px - gbCenterX;
			fVy = curPnt[m].py - gbCenterY;
			fVz = curPnt[m].pz - gbCenterZ;

			gzwrite(fo, &fVx, sizeof(float));
			gzwrite(fo, &fVy, sizeof(float));
			gzwrite(fo, &fVz, sizeof(float));
		}
	}

	cObjType = 4; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = 0;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	
	cObjType = 2; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = totalNormal * 3;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		char bNx, bNy, bNz;

		for(m=0;m<totalNormal;m++)
		{
			bNx = curPnt[m].nx; 
			bNy = curPnt[m].ny; 
			bNz = curPnt[m].nz;
			
			gzwrite(fo, &bNx, sizeof(char));
			gzwrite(fo, &bNy, sizeof(char));
			gzwrite(fo, &bNz, sizeof(char));
		}
	}

	cObjType = 3; iNumOfObjProps = 0; iNumOfObjElems = 1; 
	iLength = (int)curText.size() * 8;
	gzwrite(fo, &cObjType, sizeof(char));
	gzwrite(fo, &iNumOfObjProps, sizeof(short));
	gzwrite(fo, &iNumOfObjElems, sizeof(short));
	gzwrite(fo, &iLength, sizeof(int));
	{
		float fTx, fTy;
		
		for(m=0;m<(int)curText.size();m++)
		{
			fTx = curText[m].fU; 
			fTy = curText[m].fV;
			
			gzwrite(fo, &fTx, sizeof(float));
			gzwrite(fo, &fTy, sizeof(float));
		}
	}

	for(unsigned int i=0;i<curTri.size();i++)
	{
		int iFlag = 0;
		for(j=0; j<(int)myCultures.size(); j++)
		{
			std::string curMtype = curTri[i].cMtype;
			if((!strcmp(curTri[i].cMtype, myCultures[j].c_str())) || (curMtype.substr(0, 7) == DEFINE_AIRPORT_STRING))
			{
				iFlag = 1; break;
			}
		}
		if(!iFlag) continue;

		cObjType = curTri[i].cObjType; iNumOfObjProps = 2; iNumOfObjElems = curTri[i].iElems; 
		gzwrite(fo, &cObjType, sizeof(char));
		gzwrite(fo, &iNumOfObjProps, sizeof(short));
		gzwrite(fo, &iNumOfObjElems, sizeof(short));

		cPropType = 0; iPropLength = strlen(curTri[i].cMtype);
		gzwrite(fo, &cPropType, sizeof(char));
		gzwrite(fo, &iPropLength, sizeof(int));
		for(k=0;k<iPropLength;k++)  gzwrite(fo, &curTri[i].cMtype[k], sizeof(char));
		cPropType = 1; iPropLength = 1;  indexType = cPropVal = curTri[i].cPropVal;
		gzwrite(fo, &cPropType, sizeof(char));
		gzwrite(fo, &iPropLength, sizeof(int));
		gzwrite(fo, &cPropVal, sizeof(char));

		for(unsigned int j=0;j<curTri[i].mElement.size();j++)
		{
			unsigned short sPointV, sPointN, sPointT;

			if((indexType == 3)||(indexType == 5)||(indexType == 9))
			{
				iLength = curTri[i].mElement[j].mVertex.size() * 4;
				gzwrite(fo, &iLength, sizeof(int));

				for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size();k++)
				{
					sPointV = curTri[i].mElement[j].mVertex[k].P;
					if(indexType == 3) sPointT = curTri[i].mElement[j].mVertex[k].N;
					if(indexType == 9) sPointT = curTri[i].mElement[j].mVertex[k].T;

					gzwrite(fo, &sPointV, sizeof(short));
					gzwrite(fo, &sPointT, sizeof(short));
				}//k
			}

			if((indexType == 7)||(indexType == 11)||(indexType == 13))
			{
				iLength = curTri[i].mElement[j].mVertex.size() * 6;
				gzwrite(fo, &iLength, sizeof(int));

				for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size();k++)
				{
					sPointV = curTri[i].mElement[j].mVertex[k].P;
					sPointN = curTri[i].mElement[j].mVertex[k].N;
					sPointT = curTri[i].mElement[j].mVertex[k].T;

					gzwrite(fo, &sPointV, sizeof(short));
					gzwrite(fo, &sPointN, sizeof(short));
					gzwrite(fo, &sPointT, sizeof(short));
				}//k
			}
		}// j
	}// i

	gzclose(fo);
	return 1;
}

/* ����BTG�ļ�����BTG�ļ����ΪTXT�ļ� */
void AnalysisBTG(char *btgName)
{
	int i, j, k, m, iCreationTime, indexType;
	unsigned short iVersion, iMagicNum, iNumOfToplvlObject;
	double gbCenterX, gbCenterY, gbCenterZ;
	char simName[24];
	FILE *fout;
	gzFile fp;

	//������ļ�
	for(i=0; i<(int)strlen(btgName); i++)	
	{
		if(btgName[i] != '.') {simName[i] = btgName[i]; } else	{break;}
	}
	simName[i] = 0;
	
	sprintf(debugFile, "%s\\%s.txt",  tempDirectory, simName);
	if((fout = fopen(debugFile, "wt"))==NULL)	return;
	if((fp = gzopen(btgName, "rb")) == NULL) 	return;
	{
		//Head
		gzread(fp, &iVersion, sizeof(short));
		gzread(fp, &iMagicNum, sizeof(short));
		gzread(fp, &iCreationTime, sizeof(int));
		gzread(fp, &iNumOfToplvlObject, sizeof(short));

		fprintf(fout, "Version: %d\nMagicNumber: 0x%4X\nCreationTime: 0x%8X\nNumber of TopLevel Object: %d\n\n", iVersion, iMagicNum, iCreationTime, iNumOfToplvlObject);

		for(m=0;m<iNumOfToplvlObject;m++)
		{
			//Object Head
			char cObjType;
			unsigned short iNumOfObjProps, iNumOfObjElems;
			int iLength;

			gzread(fp, &cObjType, sizeof(char));
			gzread(fp, &iNumOfObjProps, sizeof(short));
			gzread(fp, &iNumOfObjElems, sizeof(short));

			fprintf(fout, "Object No.: %d\nObject Type: %d\nNumber of Object Properties : %d\nNumber of Object Element : %d\n", m, cObjType, iNumOfObjProps, iNumOfObjElems);

			if(iNumOfObjProps)
			{
				char cPropType;
				char sProps[64];
				int iPropLength;

				memset(sProps, 0, 64);
				for(i=0;i<iNumOfObjProps;i++)
				{
					gzread(fp, &cPropType, sizeof(char));
					gzread(fp, &iPropLength, sizeof(int));
					switch(cPropType)
					{
					case 0:		//Material
						{
							gzread(fp, sProps, sizeof(char)*iPropLength);
							fprintf(fout, "PropNo. %d  Material  : %s\n", i, sProps);
						}
						break;
					case 1:		//Index Types
						{
							gzread(fp, sProps, sizeof(char)*iPropLength);
							fprintf(fout, "PropNo. %d  IndexTypes: ", i);
							for(j=0;j<iPropLength;j++)
							{
								fprintf(fout, "%d, ", sProps[j]);
							}
							fprintf(fout, "\n");
							indexType = sProps[0];
						}
						break;
					}
				}
			}

			for(k=0;k<iNumOfObjElems;k++)
			{
				gzread(fp, &iLength, sizeof(int));
				fprintf(fout, "Length: %d\n", iLength);

				switch(cObjType)
				{
				case 0:		// Bounding Sphere  
					{
						double dXpart, dYpart, dZpart;
						float  fRadiusofBS;

						gzread(fp, &dXpart, sizeof(double));
						gzread(fp, &dYpart, sizeof(double));
						gzread(fp, &dZpart, sizeof(double));
						gzread(fp, &fRadiusofBS, sizeof(float));

						fprintf(fout, "Xpart: %g Ypart: %g Zpart %g RadiusOfBS : %f\n", dXpart, dYpart, dZpart, fRadiusofBS);

						gbCenterX = dXpart; gbCenterY = dYpart; gbCenterZ = dZpart;

					}
					break;
				case 1:		// Vertex List  
					{
						int iNumberOfVertex = iLength / 12;
						float fVx, fVy, fVz;
						Carton cart;
						Geodesy geod;

						j=0;
						for(i=0;i<iNumberOfVertex;i++)
						{
							gzread(fp, &fVx, sizeof(float));
							gzread(fp, &fVy, sizeof(float));
							gzread(fp, &fVz, sizeof(float));

							cart.x = fVx + gbCenterX; cart.y = fVy + gbCenterY; cart.z = fVz + gbCenterZ;
							CartToGeod(&cart, &geod);

							fprintf(fout, "{%13.8f, %13.8f, %13.8f}  ==  {%15.6f, %15.6f, %15.6f}--- %4d\n ", geod._lon, geod._lat, geod._elevation, cart.x, cart.y, cart.z, i); j++;
							if(j==10) {fprintf(fout, "\n"); j=0;}
						}
					}	
					fprintf(fout, "\n");
					break;
				case 2:		// Normal List  
					{
						int iNumberOfNormal = iLength / 3;
						char bNx, bNy, bNz;

						j=0;
						for(i=0;i<iNumberOfNormal;i++)
						{
							gzread(fp, &bNx, sizeof(char));
							gzread(fp, &bNy, sizeof(char));
							gzread(fp, &bNz, sizeof(char));
							fprintf(fout, "{%4d, %4d, %4d}, ", bNx, bNy, bNz); j++;
							if(j==10) {fprintf(fout, "\n"); j=0;}
						}
					}	
					fprintf(fout, "\n");
					break;
				case 3:		// Texture Coordinates List  
					{
						int iNumberOfTexture = iLength / 8;
						float fTx, fTy;

						j=0;
						for(i=0;i<iNumberOfTexture;i++)
						{
							gzread(fp, &fTx, sizeof(float));
							gzread(fp, &fTy, sizeof(float));
							fprintf(fout, "{%13.6f, %13.6f}, ", fTx, fTy); j++;
							if(j==10) {fprintf(fout, "\n"); j=0;}
						}
					}	
					fprintf(fout, "\n");
					break;
				case 4:		// Color List  
					{
						int iNumberOfColor = iLength / 16;
						float fVR, fVG, fVB, fVA;

						j=0;
						for(i=0;i<iNumberOfColor;i++)
						{
							gzread(fp, &fVR, sizeof(float));
							gzread(fp, &fVG, sizeof(float));
							gzread(fp, &fVB, sizeof(float));
							gzread(fp, &fVA, sizeof(float));
							fprintf(fout, "{%13.6f, %13.6f, %13.6f, %13.6f}, ", fVR, fVG, fVB, fVA); j++;
							if(j==10) {fprintf(fout, "\n"); j=0;}
						}
					}	
					fprintf(fout, "\n");
					break;
				case 9:		// Points  d
				case 10:	// Individual Triangles  
				case 11:	// Triangle Strips  
				case 12:	// Triangle Fans  
					{
						int iNumberOfPoints;
						unsigned short sPointV, sPointN, sPointT;

						if((indexType == 3)||(indexType == 5)||(indexType == 9))
						{
							iNumberOfPoints = iLength / 4;

							j=0;
							for(i=0;i<iNumberOfPoints;i++)
							{
								gzread(fp, &sPointV, sizeof(short));
								gzread(fp, &sPointN, sizeof(short));
								fprintf(fout, "{%3d,%3d} ", sPointV, sPointN); j++;
								if(j==10) {fprintf(fout, "\n"); j=0;}
							}
						}

						if((indexType == 7)||(indexType == 11)||(indexType == 13))
						{
							iNumberOfPoints = iLength / 6;

							j=0;
							for(i=0;i<iNumberOfPoints;i++)
							{
								gzread(fp, &sPointV, sizeof(short));
								gzread(fp, &sPointN, sizeof(short));
								gzread(fp, &sPointT, sizeof(short));

								fprintf(fout, "{%3d,%3d,%3d} ", sPointV, sPointN, sPointT); j++;
								if(j==10) {fprintf(fout, "\n"); j=0;}
							}
						}
						fprintf(fout, "\n");
					}	
					break;
				}
			}
			fprintf(fout, "\n");
		}
	}
	gzclose(fp);
	fclose(fout);
}

/* ת�����ж�������꣬����ͨ�û���������ר�õĹ��� */
int ChangeVertexInDataStru(Geodesy *source)
{
	Geodesy model, geod, dtgeo;
	Carton  src, mdl, delta, cart;
	
	/* ��������ֱ������ϵ֮��Ĳ�ֵ */
	model._lon = 108.751861; model._lat = 34.446863; model._elevation = 0.0;
	GeodToCart(&model, &mdl); GeodToCart(source, &src);
	delta.x = src.x - mdl.x; delta.y = src.y - mdl.y; delta.z = src.z - mdl.z; 
	dtgeo._lon = source->_lon - model._lon; dtgeo._lat = source->_lat - model._lat;  dtgeo._elevation = source->_elevation - model._elevation; 
	
	/* ���ݲ�ֵ���������� */
	gbCenterX += delta.x; gbCenterY += delta.y; gbCenterZ += delta.z; 
	for(int m=0;m<(int)curPnt.size();m++)
	{
		curPnt[m].lon += dtgeo._lon; 	curPnt[m].lat += dtgeo._lat; 	curPnt[m].elev += dtgeo._elevation; 
		geod._lon = curPnt[m].lon;		geod._lat = curPnt[m].lat;		geod._elevation = curPnt[m].elev;
		GeodToCart(&geod, &cart);
		curPnt[m].px = cart.x; 			curPnt[m].py = cart.y; 			curPnt[m].pz = cart.z;
	}
	
	return 1;
}

/**********************************************************
 *
 *  ���¹�����������ͨ�û�����ʱ����á�
 *
 **********************************************************/

/* �ҳ�ȹ�ߵ㼯�� */
std::vector<int> Skirts;						//ȹ�ߵ���ż���
std::vector<SimEdge> skirtEdge, sEbk;			//ȹ�ߵıߵļ���

/* ����ȹ�ߵĶ�����ţ��洢����Skirts���档 */
int FindSkirtVertex(float elev)
{
	for(int m=0;m<(int)curPnt.size();m++)
	{
		if(fabs(curPnt[m].elev - elev) > 10.0)
			Skirts.push_back(m);
	}
	return 1;
}

/* ���һ�����Ƿ���ȹ�ߵ㼯���ڣ����򷵻�0�������򷵻�1�� */
int CheckElementValid(int item)
{
	int iFlag = 1;
	for( std::vector<int>::iterator p= Skirts.begin(); p != Skirts.end(); p++ )
	{
		if(*p == item) { iFlag = 0; break; }
	}
	return iFlag;
}

/* ���һ���������ǲ���ȹ�������Σ���������ε������㶼��ȹ�ߵĵ㣬����һ����ȹ�������Ρ�
   �����ǰ��������ȹ�������Σ���ȹ�ߵ��Ǹ��ߴ��뵽ȹ�߼���skirtEdge�С�sEbk�Ǳ����õġ� */
int CheckTriangleToFindSkirtEdge(int a, int b, int c)
{
	int x, y, z;
	
	if(CheckElementValid(a) && CheckElementValid(b) && CheckElementValid(c))
	{
	}else{

		/* ��������Σ��Ƿ�����������ȹ�ߵĵ� */
		if( (!CheckElementValid(a)) && (!CheckElementValid(b)) )
		{
			x = a; y = b; z = c;
		}
		else if( (!CheckElementValid(a)) && (!CheckElementValid(c)) )
		{
			x = a; y = c; z = b;
		}
		else if( (!CheckElementValid(b)) && (!CheckElementValid(c)) )
		{
			x = b; y = c; z = a;
		}else{
			return 0;
		}

		/* ͳһ��� */
		osg::Vec3d vx; vx.set(curPnt[x].px, curPnt[x].py, curPnt[x].pz);		
		osg::Vec3d vy; vy.set(curPnt[y].px, curPnt[y].py, curPnt[y].pz);	
		if(vx != vy)
		{
			SimEdge se;	 se.A = vx; se.B = vy; se.x = x; se.y = y; se.z = z;
			skirtEdge.push_back(se);
			sEbk.push_back(se);
		}

	}
	return 1;
}

/* ��DataStruȥ��ȹ�ߵ��� */
int RemoveSkirtsFromDataStru()
{
	short a, b, c;
	
	for(unsigned int i=0; i<curTri.size(); i++)
	{
		/* type == 9Ϊ����� */
		if(curTri[i].cObjType == 9) continue;

		/* type == 10Ϊ���������� */
		if(curTri[i].cObjType == 10)
		{ 
			for(unsigned int j=0; j<curTri[i].mElement.size(); j++)
			{
				for(unsigned int m=0; m<curTri[i].mElement[j].mVertex.size(); m+=3)
				{
					a = curTri[i].mElement[j].mVertex[m + 0].P;
					b = curTri[i].mElement[j].mVertex[m + 1].P;
					c = curTri[i].mElement[j].mVertex[m + 2].P;
					CheckTriangleToFindSkirtEdge(a, b, c);
				}
			}
		}
		
		/* type == 11Ϊ����stripe */
		if(curTri[i].cObjType == 11)
		{ 
			for(unsigned int j=0; j<curTri[i].mElement.size(); j++)
			{
				a = curTri[i].mElement[j].mVertex[0].P;
				b = curTri[i].mElement[j].mVertex[1].P;
				for(unsigned int m=2; m<curTri[i].mElement[j].mVertex.size(); m+=1)
				{
					c = curTri[i].mElement[j].mVertex[m].P;
					CheckTriangleToFindSkirtEdge(a, b, c);
					a = b; b = c;
				}
			}
		}

		/* type == 12Ϊ����pan */
		if(curTri[i].cObjType == 12)
		{ 
			for(unsigned int j=0; j<curTri[i].mElement.size(); j++)
			{
				a = curTri[i].mElement[j].mVertex[0].P;
				b = curTri[i].mElement[j].mVertex[1].P;
				for(unsigned int m=2; m<curTri[i].mElement[j].mVertex.size(); m+=1)
				{
					c = curTri[i].mElement[j].mVertex[m].P;
					CheckTriangleToFindSkirtEdge(a, b, c);
					b = c;
				}
			}
		}
	}

	return 1;
}

/* ֱ�Ӹ���ȹ�ߵı߼�����Χ���õ�һ���߻�LINE_LOOP�����صĵ㼯�ϴ���coords���� */
int CalculatePerimeter( osg::Vec3dArray &coords )
{
	/* ����ȹ�߼��������ҳ�ֻ����һ�εıߣ���һ�����ڶ���������ɵ������棬ֻ����һ�εıߣ��϶�����Χ�� */
    SimEdge e = skirtEdge.front();
    osg::Vec3d A = e.A;
    osg::Vec3d B = e.B;
    do {
        for( std::vector<SimEdge>::iterator p = skirtEdge.begin(); p != skirtEdge.end(); p++ )
        {
            if( p->has(B) && !(e == *p) )
            {
                coords.push_back( A );
                A = B;
                if( B == p->A )
                    B = p->B;
                else
                    B = p->A;

                e = *p;
                break;
            }
        }
    }while( !(e == skirtEdge.front()) );
    coords.push_back(A);
    return 1;
}

/* �ռ�ȫ����LINE_LOOP */
osg::ref_ptr<osg::Geode> GetAllLineLoopFromSkirtEdge()
{
	/* ����һ��Geode */
	osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
	gnode->setName("River");	
	
	while(skirtEdge.size() > 0)
	{
		/* ���㵱ǰ��Χ*/
		osg::ref_ptr<osg::Vec3dArray> perimeter = new osg::Vec3dArray;
		CalculatePerimeter( *perimeter.get() );
		
		/* ����һ��Geometry */
		osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
		gnode->addDrawable(gm);
		gm->setVertexArray(perimeter);
		gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, perimeter->getNumElements()));
		
		/* �Ѽ��������Χ�ߴ�skirtEdge�޳���ȥ */
		for( osg::Vec3dArray::iterator v = perimeter->begin(); v != perimeter->end(); v++ )
		{
			for( std::vector<SimEdge>::iterator e = skirtEdge.begin(); e != skirtEdge.end(); e++ )
			{
				if(e->has(*v))
				{
					skirtEdge.erase(e);
					break;
				}
			}

			/* �ҵ�ȹ�ߵ��Ӧ�����ϵĵ� */
			for( std::vector<SimEdge>::iterator e = sEbk.begin(); e != sEbk.end(); e++ )
			{
				int z = e->get(*v);
				if((e->has(*v))&&(z != -1))
				{
					v->set(curPnt[z].px, curPnt[z].py, curPnt[z].pz);
					break;
				}
			}

		}//for
	}//while

	/* �ͷ�������ڴ� */
	skirtEdge.clear();
	sEbk.clear();
	Skirts.clear();
	return gnode.release();
}

/**********************************************************
 *
 *  OSG ��һЩ���̣� ���߻�����߽�
 *
 **********************************************************/

/* ���TerrainProcessingUtils���̡�*/
int computePerimeter( std::vector<SimEdge> &edges, osg::Vec3Array &coords )
{
    for( std::vector<SimEdge>::iterator p = edges.begin(); p != edges.end(); p++ )
    {
        std::vector<SimEdge>::iterator q = p;
        q++;
        while(q != edges.end() )
        {
            if( *p == *q )
            {
                p->keep = false;
                q->keep = false;
            }
            q++;
        }
    }

    std::vector<SimEdge> kept;
    for( std::vector<SimEdge>::iterator p = edges.begin(); p != edges.end(); p++ )
    {
        if( p->keep )
            kept.push_back( *p );
    }
    if(kept.size() == 0) return 0;

    SimEdge e = kept.front();
    osg::Vec3 A = e.A;
    osg::Vec3 B = e.B;
    do {
        for( std::vector<SimEdge>::iterator p = kept.begin(); p != kept.end(); p++ )
        {
            if( p->has(B) && !(e == *p) )
            {
                coords.push_back( A );
                A = B;
                if( B == p->A )
                    B = p->B;
                else
                    B = p->A;

                e = *p;
                break;
            }
        }
    }while( !(e == kept.front()) );
    coords.push_back(A);
    return 1;
}

/* ��֪һ����άƽ�漸���壬��������������Χ�ߵ��߻�LineLoop�� */
osg::ref_ptr<osg::Geode> CalcuPerimeterLINELOOP(osg::ref_ptr<osg::Geode> gnode)
{
	/* ����gnode��ÿһ��geometry����Ҫ��LINE_LOOP */
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	
    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    {
        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
        if( geometry.valid() )
        {
		    std::vector<SimEdge>   edges;
		    osg::TriangleFunctor<SimMakeEdgeList> mel;
		    mel.setEdgeVector(&edges);
		    geometry->accept( mel );

			while(edges.size() > 0)
			{
	        	/* �����һ��LINE_LOOP */
				osg::ref_ptr<osg::Vec3Array> perimeter = new osg::Vec3Array;
				if(computePerimeter( edges, *perimeter.get()) == 0) break;

				/* ��perimeter��LINE_LOOP���뵽һ��geometry���� */
				if (perimeter->size()>=3)
				{
					osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
					gd->addDrawable(gm);
					gm->setVertexArray( perimeter.get() );
					gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, perimeter->size()) );
				}

				/* ������keep����ԭΪtrue */
			    for( std::vector<SimEdge>::iterator p = edges.begin(); p != edges.end(); p++ ) 
					p->keep = true;

				
				/* ��edges�����߼��������perimeters */
				for( osg::Vec3Array::iterator v = perimeter->begin(); v != perimeter->end(); v++ )
				{
					for( std::vector<SimEdge>::iterator e = edges.begin(); e != edges.end(); e++ )
					{
						if(e->has(*v)) e->keep = false;
					}
				}

				/* ��keep == false�Ķ�ȥ�� */
			    std::vector<SimEdge> kept;
			    for( std::vector<SimEdge>::iterator p = edges.begin(); p != edges.end(); p++ )
			    {
			        if( p->keep )
			            kept.push_back( *p );
			    }
			    edges.clear();				    
				for( std::vector<SimEdge>::iterator e = kept.begin(); e != kept.end(); e++ )
				{
					edges.push_back( *e );
				}
				
			}	//while
        }	//if
    }	//for
    
    
    /* ����gd��ÿһ��LINE_LOOP����Ҫ�ӽ��߼��Ͽ����Ƿ����ظ��ıߡ���������ܵ�û���رߵ�Line_Loop,*/
    std::vector<SimEdge> alledges;	alledges.clear();
    for( unsigned int n = 0; n < gd->getNumDrawables(); n++ )
    {
        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gd->getDrawable(n));
		osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray() );
        if( geometry.valid() )
        {
		    for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
		    {
		        osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
		        
				switch( pset->getType() )
				{
				    case osg::PrimitiveSet::DrawArraysPrimitiveType:
				    {
				        osg::DrawArrays *da = dynamic_cast<osg::DrawArrays *>(pset.get());
				        int first = da->getFirst();
				        unsigned int count = da->getCount();
				
				        switch(pset->getMode() )
				        {
				            case osg::PrimitiveSet::LINE_LOOP:
				            {
								for( unsigned int ii = 1; ii < count; ii++ )
								{
									osg::Vec3 vx, vy;
									vx = (*coords)[first + ii - 1];
									vy = (*coords)[first + ii - 0];
									SimEdge se;	 se.A = vx; se.B = vy; se.x = 0; se.y = 0; se.z = 0;
									alledges.push_back(se);
								}
				            }
				            break;
				        }
				    }
				    break;
				}
		    }
		}
	}	    
	
	gd->removeDrawables(0, gd->getNumDrawables());
	while(alledges.size() > 0)
	{
    	/* �����һ��LINE_LOOP */
		osg::ref_ptr<osg::Vec3Array> perimeter = new osg::Vec3Array;
		if(computePerimeter( alledges, *perimeter.get()) == 0) break;

		/* ��perimeter��LINE_LOOP���뵽һ��geometry���� */
		if (perimeter->size()>=3)
		{
			osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
			gd->addDrawable(gm);
			gm->setVertexArray( perimeter.get() );
			gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, perimeter->size()) );
		}

		/* ������keep����ԭΪtrue */
	    for( std::vector<SimEdge>::iterator p = alledges.begin(); p != alledges.end(); p++ ) 
			p->keep = true;
		
		/* ��alledges�����߼��������perimeters */
		for( osg::Vec3Array::iterator v = perimeter->begin(); v != perimeter->end(); v++ )
		{
			for( std::vector<SimEdge>::iterator e = alledges.begin(); e != alledges.end(); e++ )
			{
				if(e->has(*v)) e->keep = false;
			}
		}

		/* ��keep == false�Ķ�ȥ�� */
	    std::vector<SimEdge> kept;
	    for( std::vector<SimEdge>::iterator p = alledges.begin(); p != alledges.end(); p++ )
	    {
	        if( p->keep )
	            kept.push_back( *p );
	    }
	    alledges.clear();				    
		for( std::vector<SimEdge>::iterator e = kept.begin(); e != kept.end(); e++ )
		{
			alledges.push_back( *e );
		}
		
	}	//while

	return gd.release();	
}


/* �Ѵ�btg�ļ�����Ľṹ��DataStru���ΪGroup */
osg::ref_ptr<osg::Group> CreateGroupFromDataStru()
{
	int offset;
	float fNx, fNy, fNz;
	
	//ȱʡһ������
	fNx = 92/127.5; fNy = -30/127.5; fNz = -57/127.5; offset = 0;

	osg::ref_ptr<osg::Group> group = new osg::Group;
	group->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);

	//����ÿ�������棬����һ�������壬��Ϊһ��drawable���뵽gnode������
	for(unsigned int i=0; i<curTri.size(); i++)
	{
		if(curTri[i].cObjType == 9) continue;

		osg::ref_ptr<osg::Geode> gnode = new osg::Geode;
		for(unsigned int j=0; j<curTri[i].mElement.size(); j++)
		{
			osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
			gnode->addDrawable(gm);

			//ѹ�붥��
			osg::ref_ptr<osg::Vec3Array> vertex = new osg::Vec3Array;
			for(unsigned int m=0; m<curTri[i].mElement[j].mVertex.size(); m++)
			{
				osg::Vec3d v; 
				v.set(curPnt[curTri[i].mElement[j].mVertex[m].P].px, curPnt[curTri[i].mElement[j].mVertex[m].P].py, curPnt[curTri[i].mElement[j].mVertex[m].P].pz);
				vertex->push_back(v);				
			}
			gm->setVertexArray(vertex);
			
			//ѹ�뷨��
			osg::ref_ptr<osg::Vec3Array> normal = new osg::Vec3Array;
			normal->push_back(osg::Vec3(fNx, fNy, fNz));
			gm->setNormalArray(normal);
			gm->setNormalBinding(osg::Geometry::BIND_OVERALL);

			//ѹ����������
			osg::ref_ptr<osg::Vec2Array> coord = new osg::Vec2Array;
			for(unsigned int m=0; m<curTri[i].mElement[j].mVertex.size(); m++)
			{
				osg::Vec2 c; c.set(curText[curTri[i].mElement[j].mVertex[m].P].fU, curText[curTri[i].mElement[j].mVertex[m].P].fV);
				coord->push_back(c);
			}
			gm->setTexCoordArray(0, coord);
			
			//ѹ����
			switch(curTri[i].cObjType)
			{
			case 10:
		    	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLES,0,curTri[i].mElement[j].mVertex.size()));
		    	break;
		    case 11:
		    	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLE_STRIP,0,curTri[i].mElement[j].mVertex.size()));
		    	break;
		    case 12:
		    	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLE_FAN,0,curTri[i].mElement[j].mVertex.size()));
		    	break;
		    }			
		}
		
		/* ������������ʲô������ʽ��������������ʽ������������߻�����ʽ�����*/
		if(1)	//(!((!strcmp(curTri[i].cMtype, "Stream")) || (!strcmp(curTri[i].cMtype, "Road"))))
		{
			gnode->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
			gnode->setName(curTri[i].cMtype);
			group->addChild(gnode.get());
			gnode.release();
		}else{
			/* �����ܵ�perimeter��LINE_LOOP */
			osg::ref_ptr<osg::Geode> gd = CalcuPerimeterLINELOOP(gnode);
			//osgDB::writeNodeFile(*gd, "E:\\GENZ.ive");
			group->addChild(gd.get());
		}
	}

	return group.release();
} 

/* У��ͨ�û���ģ�塣 */
int VerifyTheModels(char *fn, char *out)
{
	int j, m;
	char meter[32];

	/* �ȶ���Btg�ļ� */
	ReadBTGToDataStruct(fn);
	
	/* ���Vertex�߶��Ƿ�Ϊ0 */
	if(fabs(curPnt[0].elev) > 10.0)
	{
		Geodesy geod;
		Carton  cart;

		for(m=0;m<(int)curPnt.size();m++)
		{
			curPnt[m].elev = 0.0; 
			geod._lon = curPnt[m].lon;		geod._lat = curPnt[m].lat;		geod._elevation = curPnt[m].elev;
			GeodToCart(&geod, &cart);
			curPnt[m].px = cart.x; 			curPnt[m].py = cart.y; 			curPnt[m].pz = cart.z;
		}
	}
	
	/* �����ƹ�Ĳ��� */
	for(unsigned int i=0; i<curTri.size(); i++)
	{
		if(curTri[i].cObjType == 9)
		{
			if((curTri[i].cMtype[0] == '0') || (curTri[i].cMtype[0] == '1') || (curTri[i].cMtype[0] == 'x') || (curTri[i].cMtype[0] == 'X'))
			{
				int len = strlen(curTri[i].cMtype);
				memset(meter, 0, sizeof(meter)); m = 0;
				for(j=4; j<len; j++) meter[m++] = curTri[i].cMtype[j];
				memset(curTri[i].cMtype, 0, sizeof(curTri[i].cMtype));
				strcpy(curTri[i].cMtype, meter);
			}
		}
	}
		
	/* д�ص�btg�ļ����� */
	WriteDataStructToBTG(out);	

	return 1;
}

/*******************************************************
 *
 *  Ѱ��һ��group��perimeter
 *  ���¹�����ͨ���û�flt�����������и�ؾ���ʱ���ã���������flt����������
 *
 *******************************************************/

std::vector<FoundNode>  allFn;			//�������в鵽�Ķ���
std::vector<string> 	gname;			//�������в鵽������
FoundNode mNode;						//���浱ǰ��flt������ֵ��

/* �������е�Geode���ҳ�ÿ��Geode�����ж�������ֵ����Сֵ��������allFn���� */
class FindMinMaxVisitor : public osg::NodeVisitor
{
public:
	FindMinMaxVisitor() : osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { _elevationFlag = false; }; 

	void apply(osg::Group &gp)
	{
		for(unsigned int i = 0; i < gp.getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp.getChild(i));
			if(gnode != NULL)
				apply(*gnode);
		}
		traverse(gp);
	}

	void apply(osg::Geode &gnode)
	{
		int idx = 0, iFlag = 0;

		osg::Geode::DrawableList dl = gnode.getDrawableList();
		if(dl.size() == 0) return;
		for(osg::Geode::DrawableList::iterator iter = dl.begin(); iter != dl.end(); iter++)
		{
			float minX, minY, maxX, maxY, minZ, maxZ;
			minX = minY = minZ = 99999999.9;  maxX = maxY = maxZ =-99999999.9;
			
			osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry*> ((*iter).get());
			osg::Vec3Array *vx = dynamic_cast<osg::Vec3Array*>(gm->getVertexArray());
			if(vx == NULL ) continue;
			for(osg::Vec3Array::iterator v = vx->begin(); v != vx->end(); v++)
			{
				if(v->x() > maxX) maxX = v->x();
				if(v->x() < minX) minX = v->x();
				if(v->y() > maxY) maxY = v->y();
				if(v->y() < minY) minY = v->y();
				if(v->z() > maxZ) maxZ = v->z();
				if(v->z() < minZ) minZ = v->z();
				if((v->z() > 10.0)&&(_elevationFlag)) { iFlag = 1; break;}
			}

			if(!iFlag)
			{
				FoundNode fn;
				fn.name  = gnode.getName();
				fn.gmIdx = idx;
				fn.maxX  = maxX;
				fn.maxY  = maxY;
				fn.maxZ  = maxZ;
				fn.minX  = minX;
				fn.minY  = minY; 
				fn.minZ  = minZ;
				allFn.push_back(fn);
			}
			idx++;
		}
		traverse(gnode);
	}

	void setElevationFlag(bool ef)
	{
		_elevationFlag = ef;
	}

private:
	bool _elevationFlag;
};

/* ��allFn���ҳ����ֵ����Сֵ�����ҳ���Ӧ��Geode����, ɸѡ�����ظ������ƣ�������gname���� */
void FindGnameInAllFn()
{
	FoundNode maX, maY, miX, miY, bFn;
	int iFlag = 0;
	
	if(allFn.size() == 0)	return;

	bFn = allFn.front();
	maX.maxX = bFn.maxX; maX.name = bFn.name;
	maY.maxY = bFn.maxY; maY.name = bFn.name;
	miX.minX = bFn.minX; miX.name = bFn.name;
	miY.minY = bFn.minY; miY.name = bFn.name;

	for(std::vector<FoundNode>::iterator p = allFn.begin(); p != allFn.end(); p++)
	{
		if(p->maxX > maX.maxX) { maX.maxX = p->maxX; maX.name = p->name;}
		if(p->maxY > maY.maxY) { maY.maxY = p->maxY; maY.name = p->name;}
		if(p->minX < miX.minX) { miX.minX = p->minX; miX.name = p->name;}
		if(p->minY < miY.minY) { miY.minY = p->minY; miY.name = p->name;}
	}

	mNode.maxX = maX.maxX;
	mNode.maxY = maY.maxY;
	mNode.minX = miX.minX;
	mNode.minY = miY.minY;
	
	/* ɸѡ */
	string firstname = maX.name;	gname.push_back(firstname);

	string othername = maY.name; iFlag = 0;
	for(std::vector<string>::iterator p = gname.begin(); p != gname.end(); p++)
	{
		if(*p == othername)	{iFlag = 1; break; }
	}
	if(!iFlag) gname.push_back(othername);

	othername = miX.name; iFlag = 0;
	for(std::vector<string>::iterator p = gname.begin(); p != gname.end(); p++)
	{
		if(*p == othername)	{iFlag = 1; break; }
	}
	if(!iFlag) gname.push_back(othername);

	othername = miY.name; iFlag = 0;
	for(std::vector<string>::iterator p = gname.begin(); p != gname.end(); p++)
	{
		if(*p == othername)	{iFlag = 1; break; }
	}
	if(!iFlag) gname.push_back(othername);

	/* �ٱ���һ��allFn, ��ӵ�м�ֵ����ȴû����gname��������Geode������Ҳ�н���gname�б��� 2012-1-12 */
	for(std::vector<FoundNode>::iterator p = allFn.begin(); p != allFn.end(); p++)
	{
		if((p->maxX == maX.maxX) || (p->maxY == maY.maxY) || (p->minX == miX.minX) || (p->minY == miY.minY))
		{
			iFlag = 0;
			for(std::vector<string>::iterator q = gname.begin(); q != gname.end(); q++)
			{	if(*q == p->name)	{iFlag = 1; break; }	}
			if(!iFlag) gname.push_back(p->name);
		}
	}
}

/* ��allFn���ҳ����ֵ����Сֵ��������mNodeȫ�ֱ������� */
void FindMaxMinInAllFn()
{
	FoundNode maX, maY, maZ, miX, miY, miZ, bFn;
	int iFlag = 0;
	
	bFn = allFn.front();
	maX.maxX = bFn.maxX;
	maY.maxY = bFn.maxY;
	maZ.maxZ = bFn.maxZ;
	miX.minX = bFn.minX;
	miY.minY = bFn.minY;
	miZ.minZ = bFn.minZ;

	for(std::vector<FoundNode>::iterator p = allFn.begin(); p != allFn.end(); p++)
	{
		if(p->maxX > maX.maxX) { maX.maxX = p->maxX;}
		if(p->maxY > maY.maxY) { maY.maxY = p->maxY;}
		if(p->maxZ > maZ.maxZ) { maZ.maxZ = p->maxZ;}
		if(p->minX < miX.minX) { miX.minX = p->minX;}
		if(p->minY < miY.minY) { miY.minY = p->minY;}
		if(p->minZ < miZ.minZ) { miZ.minZ = p->minZ;}
	}
	
	mNode.maxX = maX.maxX;
	mNode.maxY = maY.maxY;
	mNode.maxZ = maZ.maxZ;
	mNode.minX = miX.minX;
	mNode.minY = miY.minY;
	mNode.minZ = miZ.minZ;
}


/* �������е�Geode���ҳ��������ֵ�Geode���棬���浽һ��Group�� */
class FindNameVisitor : public osg::NodeVisitor
{
public:
	FindNameVisitor() : osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) {}; 

	void apply(osg::Geode &gnode)
	{
		for(std::vector<string>::iterator p = gname.begin(); p != gname.end(); p++)
		{
			if((gnode.getName() == (*p)) && (group.valid()))
			{
				group->addChild(&gnode);	
			}
		}
		
		traverse(gnode);
	}

	void setGroup(osg::Group *gp)
	{
		group = gp;
	}
private:
	osg::ref_ptr<osg::Group> group;
};

#define THHOLD 1.0
/* ��ȡһ���߶εĴ�ֱƽ�����ϵ����������㣬���ж����Ǻ���֪���Ƿ��ཻ */
bool DetectIntersectOfEdge(SimEdge *se, osg::ref_ptr<osg::Node> terrain)
{
	/* �������֪�߶ε��е� */
	osg::Vec3d mid = (se->A + se->B) / 2;
	osg::Vec3d V[2];
	
	/* �������֪�߶ε�б�� */
	if(se->A[0] != se->B[0])
	{
		float k = (se->A[1] - se->B[1]) / (se->A[0] - se->B[0]);
		
		/* ������ֱ�ֱཻ�ߵ�б����˻�Ϊ-1��k1*k2=-1. */
		float kv;
		if( k != 0.0) 	kv = -(1.0 / k);
		
		/* �����ֱƽ������������ */
		V[0].x() = mid.x() + THHOLD; V[1].x() = mid.x() - THHOLD;
		if( k != 0.0)
		{
			V[0].y() = kv * (V[0].x() - mid.x()) + mid.y(); 
			V[1].y() = kv * (V[1].x() - mid.x()) + mid.y(); 
		}else{
			V[0].x() = V[1].x() = mid.x();
			V[0].y() = mid.y() + THHOLD; V[1].y() = mid.y() - THHOLD;
		}
		V[0].z() = mid.z(); V[1].z() = mid.z();
	}else{
		V[0].x() = mid.x() + THHOLD; V[1].x() = mid.x() - THHOLD;
		V[0].y() = V[1].y() = mid.y();
		V[0].z() = mid.z(); V[1].z() = mid.z();
	}

	/* �ж��ཻ */
	int iFlag = 0;
	for(unsigned int i = 0; i < 2; i++)
	{
	    osg::Vec3d p0 = V[i] + osg::Vec3d(0,0,10000.0);
	    osg::Vec3d p1 = V[i] + osg::Vec3d(0,0,-10000.0);
	
	    osg::ref_ptr<osg::LineSegment> seg = new osg::LineSegment(p0,p1);
	    if (!seg->valid())
	    {
	        i = i;
	    }
	
	    osgUtil::IntersectVisitor iv;
	    iv.addLineSegment( seg.get() );
	    terrain->accept( iv );
	
	    if( iv.hits() )
	    {
	    	iFlag++;
	    }		
	}
	
	if(iFlag==1) return false;
	if(iFlag==2) return true;
	return true;
}

/* ��������peri��������б����Ϊһ��ive�ļ������Բ鿴�м����õĺ�����*/
void OutPutEdgesIve(std::vector<SimEdge> peri, char *fn)
{
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	osg::ref_ptr<osg::Geometry> linesGeom = new osg::Geometry;
	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	for( std::vector<SimEdge>::iterator p = peri.begin(); p != peri.end(); p++ )
	{
		osg::Vec3 va = p->A;   	vx->push_back(va);
		osg::Vec3 vb = p->B;	vx->push_back(vb);
	}
	linesGeom->setVertexArray(vx);
	// set the colors as before, plus using the above
	osg::Vec4Array* colors = new osg::Vec4Array;
	colors->push_back(osg::Vec4(1.0f,1.0f,0.0f,1.0f));
	linesGeom->setColorArray(colors);
	linesGeom->setColorBinding(osg::Geometry::BIND_OVERALL);
	// set the normal in the same way color.
	osg::Vec3Array* normals = new osg::Vec3Array;
	normals->push_back(osg::Vec3(0.0f,-1.0f,0.0f));
	linesGeom->setNormalArray(normals);
	linesGeom->setNormalBinding(osg::Geometry::BIND_OVERALL);
	// This time we simply use primitive, and hardwire the number of coords to use 
	// since we know up front,
	linesGeom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,vx->getNumElements()));
	gd->addDrawable(linesGeom);
	osgDB::writeNodeFile(*gd, fn);    	
}

/* ���㼯�����Ϊһ��LINE_LOOP��ive�ļ������Բ鿴�м����õĺ��� */
void OutPutLineLoopIve(osg::ref_ptr<osg::Vec3Array> coords, char *fn)
{
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	gd->addDrawable(gm);
	gm->setVertexArray( coords.get() );
	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, coords->size()) );
	// set the colors as before, plus using the above
	osg::Vec4Array* colors = new osg::Vec4Array;
	colors->push_back(osg::Vec4(1.0f,1.0f,0.0f,1.0f));
	gm->setColorArray(colors);
	gm->setColorBinding(osg::Geometry::BIND_OVERALL);
	osgDB::writeNodeFile(*gd, fn);
}
void OutPutLineLoopIve2(osg::ref_ptr<osg::Vec3dArray> coords, char *fn)
{
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	gd->setName("River");
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	gd->addDrawable(gm);
	gm->setVertexArray( coords.get() );
	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, coords->size()) );
	// set the colors as before, plus using the above
	osg::Vec4Array* colors = new osg::Vec4Array;
	colors->push_back(osg::Vec4(1.0f,1.0f,0.0f,1.0f));
	gm->setColorArray(colors);
	gm->setColorBinding(osg::Geometry::BIND_OVERALL);
	osgDB::writeNodeFile(*gd, fn);
}


/*******************************************************
 *
 *  ���¹�����ͨ���û�flt�����������и�ؾ���ʱ���ã���������flt����������
 *
 *******************************************************/
//#define USE_TOTALEDGES
std::vector<SimEdge>  totalEdges;
bool foundAround = false;

/* �������е�Geode���ҳ�ÿ��Geode�����еıߣ����ж��Ƿ��б߽�ıߣ�����У����뵽allFn */
class FindEdgeVisitor : public osg::NodeVisitor
{
public:
	FindEdgeVisitor() : osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { gmEdges.clear(); }; 

	/* Add by Galen 2013-03-15 */
	void apply(osg::Geode &gd)
	{
		if ((gd.getName().substr(0, strlen("Around")) == "Around") && (aroundGD != NULL))
		{
			for(unsigned int j = 0; j < gd.getNumDrawables(); j++)
			{
				osg::ref_ptr<osg::Drawable> dr = dynamic_cast<osg::Drawable*> (gd.getDrawable(j));
				if(dr.valid()) aroundGD->addDrawable(dr);
			}
		}
	}
	/* Add by Galen 2013-03-15 End */

	void apply(osg::Group &gp)
	{
		if (gp.getName().substr(0, strlen("Lights")) == "Lights")
		{
			char aptLights[256];
			sprintf(aptLights, "%s\\AptLight.ive", tempDirectory);
			osgDB::writeNodeFile(gp, aptLights);
		}

		if ((gp.getName().substr(0, strlen("Around")) == "Around") && (aroundGD != NULL))
		{
			foundAround = true;
			for(unsigned int i = 0; i < gp.getNumChildren(); i++)
			{
				osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp.getChild(i));
				if(gnode != NULL)
				{
#ifdef	USE_TOTALEDGES
					doApply(*gnode);
#else		//USE_TOTALEDGES
					for(unsigned int j = 0; j < gnode->getNumDrawables(); j++)
					{
						osg::ref_ptr<osg::Drawable> dr = dynamic_cast<osg::Drawable*> (gnode->getDrawable(j));
						if(dr.valid()) aroundGD->addDrawable(dr);
					}
#endif	//USE_TOTALEDGES
				}
			}
		}
		traverse(gp);
	}

#ifdef	USE_TOTALEDGES
	bool isPureTrangle(osg::Geometry *g)
	{
	    for( unsigned int i = 0; i < g->getNumPrimitiveSets(); i++ )
	    {
	        osg::ref_ptr<osg::PrimitiveSet> pset = g->getPrimitiveSet(i);
	        if(	( pset->getMode() == osg::PrimitiveSet::QUADS) ||
	        	( pset->getMode() == osg::PrimitiveSet::QUAD_STRIP) ||
	        	( pset->getMode() == osg::PrimitiveSet::LINES) ||
	        	( pset->getMode() == osg::PrimitiveSet::LINE_STRIP) ||
	        	( pset->getMode() == osg::PrimitiveSet::LINE_LOOP) ||
	        	( pset->getMode() == osg::PrimitiveSet::POLYGON) )
	        		return false;
		}
		return true;
	}

	void doApply(osg::Geode &gnode)
	{
		int idx = 0, iFlag = 0;

		osg::Geode::DrawableList dl = gnode.getDrawableList();
		if(dl.size() == 0) return;
		for(osg::Geode::DrawableList::iterator iter = dl.begin(); iter != dl.end(); iter++)
		{
			gmEdges.clear();

			/* �������ǰGeometry�����еı� */
			osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry*> ((*iter).get());
			
			/* �����ǰ��Geometry���Ǵ�����������ɣ��������ı��εȵȣ��򲻱ؼ��㡣*/
			/* ������֤������ʹ������QUADS��Ҳ��Ӱ�������������һ����ȡ�������������ɹ��ĸ��ʡ�*/
#ifdef	PICKOUTQUADS
			int iFlag = 0;
	    for( unsigned int i = 0; i < gm->getNumPrimitiveSets(); i++ )
	    {
	        osg::ref_ptr<osg::PrimitiveSet> pset = gm->getPrimitiveSet(i);
	        if(!(	( pset->getMode() == osg::PrimitiveSet::QUADS) ||
	        		( pset->getMode() == osg::PrimitiveSet::QUAD_STRIP) ||
	        		( pset->getMode() == osg::PrimitiveSet::LINES) ||
	        		( pset->getMode() == osg::PrimitiveSet::LINE_STRIP) ||
	        		( pset->getMode() == osg::PrimitiveSet::LINE_LOOP) ||
	        		( pset->getMode() == osg::PrimitiveSet::POLYGON) )	)
	        	{		iFlag = 1; break;  }
			}
			if(!iFlag) continue;
#endif	//PICKOUTQUADS
			//if(isPureTrangle(gm) == false) continue;

	    if( gm.valid() )
	    {
			    osg::TriangleFunctor<SimMakeEdgeList> mel;
			    mel.setEdgeVector(&gmEdges);
			    gm->accept( mel );
			}
			
			if(gmEdges.size() == 0) continue;
			/* �ж�ÿ���ߣ��Ƿ��Ǳ߽�ı� */
			for( std::vector<SimEdge>::iterator p = gmEdges.begin(); p != gmEdges.end(); p++ )
			{
				SimEdge se = *p;
				
				/* �����һ�������ڱ߽��ڵıߣ��Ǹ�����������������бߣ����Ǳ߽��ڵıߡ�*/
				if(DetectIntersectOfEdge(&se, terrain) == false)
				{
					totalEdges.push_back(se);
				}
			}
		}
		traverse(gnode);
	}
#endif	//USE_TOTALEDGES

	void setTerrain(osg::Node *te)
	{
		terrain = te;
	}

	void setGeode(osg::Geode *gd)
	{
		aroundGD = gd;
	}

private:
	std::vector<SimEdge>  gmEdges;
	osg::ref_ptr<osg::Node> terrain;
	osg::ref_ptr<osg::Geode> aroundGD;
};


/* ��ȡһ������FLT�ļ������ƹ�Lights�飩��������������������ĵƹ���Ϣ��ʱ�ļ� AptLight.ive */
int ReadFltAndCreateAptLight(char *fname)
{
	osg::ref_ptr<osg::Node> te = osgDB::readNodeFile(fname);

	/* �ҳ��ƹ���ʱ�ļ���*/
	FindEdgeVisitor fev;
	fev.setTerrain(te);
	fev.setGeode(NULL);
	te->accept(fev);

	return 1;
}



/* ��ȡһ������FLT�ļ�����������Around�飩��������������������ı�ԵLineLoop����λ�û������õ���Ҫ���� */
int ReadFltAndVistIt(char *fname, osg::ref_ptr<osg::Geode> gnode)
{
	osg::ref_ptr<osg::Node> te = osgDB::readNodeFile(fname);
	
	/* ���Ȳ��Ҽ�ֵ���鵽��ֵ�ͻ�鵽�߽� */
	allFn.clear();
	FindMinMaxVisitor fmmv;
	fmmv.setElevationFlag(true);
	te->accept(fmmv);
	
	/* �ҳ������ĸ��߽ǵ�ƽ���Geode�������б��ҳ���ֵmNode���Ա��������ж� */
	gname.clear();
	FindGnameInAllFn();
	
	/* �ҳ����б�Ե�ıߣ��ô�ֱƽ�ֵ��жϷ���������֤������ʹ������QUADS��Ҳ��Ӱ���������*/
	totalEdges.clear();
	FindEdgeVisitor fev;
	fev.setTerrain(te);
	fev.setGeode(gnode);
	te->accept(fev);

#ifdef	USE_TOTALEDGES

#if 0
	sprintf(debugFile, "%s\\Flt0.ive", tempDirectory);
	OutPutEdgesIve(totalEdges, debugFile);
#endif

	if(totalEdges.size() == 0)	return 0;

	/* �����е�Edge�����ֻ�е���ͬ���Edge���ӣ�����������Edge*/
	int lastNum, currNum = 0;
	std::vector<SimEdge> peri; 
	do {
		peri.clear(); lastNum = currNum;
		
		/* �ҳ��������ӻ����ޱ����ӵ�Edge��*/
		for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
		{
			int iFlagA, iFlagB;
			iFlagA = iFlagB = 0;
	
			osg::Vec3d A = p->A;
			osg::Vec3d B = p->B;
	
			for( std::vector<SimEdge>::iterator q = totalEdges.begin(); q != totalEdges.end(); q++ )
			{
	    		if ((!( *p == *q )) && ( q->has(A)))	iFlagA = 1;
	    		if ((!( *p == *q )) && ( q->has(B)))	iFlagB = 1;
			}
			if((iFlagA) && (iFlagB))
			{
				SimEdge se = *p;
				peri.push_back(se);
			}
		}
		currNum = peri.size();
		totalEdges.clear();	for( std::vector<SimEdge>::iterator q = peri.begin(); q != peri.end(); q++ ) totalEdges.push_back(*q);
		
	}while(lastNum != currNum);

	if(peri.size() == 0)	return 0;

#if 0
	sprintf(debugFile, "%s\\Flt1.ive", tempDirectory);
	OutPutEdgesIve(peri, debugFile);
#endif

	/* �������бߵ����ж˵㣬�����һ���������߻�����ôÿ���˵�ֻ����2�Ρ�*/
	/* �����ĳ���˵�ʹ�õĴ�������2�Σ���ô˵��������Ƕ�׻���*/
	std::vector<SimEdge> kept;
	for( std::vector<SimEdge>::iterator p = peri.begin(); p != peri.end(); p++ )
	{
		if(( p->keep ) && (p->A != p->B ))
			kept.push_back( *p );
	}
	
	std::vector<SimEdge> hold;
	int lastEdgeSize = 0;
	while (1)
	{
		unsigned int i, j;

		/* ���߶��������еĵ㶼���ظ����ҳ�����*/
		osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;
		SimEdge e = kept.front();
		coords->push_back( e.A );
		coords->push_back( e.B );
		for( std::vector<SimEdge>::iterator p = kept.begin() + 1; p != kept.end(); p++ )
		{
			osg::Vec3 A = p->A;
			osg::Vec3 B = p->B;
			int iFlagA = 0, iFlagB = 0;
			for(osg::Vec3Array::iterator r = coords->begin(); r != coords->end(); r++)
			{
				if(A == *r)	iFlagA = 1;
				if(B == *r)	iFlagB = 1;
			}
			if(!iFlagA)	coords->push_back(A);
			if(!iFlagB)	coords->push_back(B);
		}

		/* ��һ�γ���Ϊ��ɾ���߻�������������Ȧ��ʹ���߻����һ�������߻���*/
		std::vector<VerIdx> OddPointsU2;
		std::vector<VerIdx> OddPointsE1;
		
		/* ��һ��Ҫ�ҳ����е���㣬���Ķ����ǣ���һ���߻��г��ֵĴ���������2�ĵ㣬һ�����߻��ϵ�ÿ����Ӧ�ö�ֻ����2�Ρ� */
		osg::ref_ptr<osg::Vec3Array> vx = coords;
		int vNum = vx->getNumElements();
		int index = 0;
		for( osg::Vec3Array::iterator r = vx->begin(); r != vx->end(); r++ )
		{
			int times = 0;
			for( std::vector<SimEdge>::iterator p = kept.begin(); p != kept.end(); p++ )
			{
	    		if(*r == p->A)	times++;
	    		if(*r == p->B)	times++;
			}
			
			if(times > 2) 
			{
				VerIdx vi; vi.Ver = *r; vi.idx = index; vi.times = times;
				OddPointsU2.push_back(vi);
			}

			if(times == 1) 
			{
				VerIdx vi; vi.Ver = *r; vi.idx = index; vi.times = times;
				OddPointsE1.push_back(vi);
			}

			index++;
		}

		/* �ڶ������������֮������µ��߶�  */
		std::vector<SimEdge> OddLineSegU2;
		if(OddPointsU2.size() > 0)
		{
			for(i=0; i<OddPointsU2.size() - 1; i++)
			{
				for(j=i+1; j<OddPointsU2.size(); j++)
				{
					SimEdge se;	 se.A = OddPointsU2[i].Ver; se.B = OddPointsU2[j].Ver; se.x = 0; se.y = 0; se.z = 0;
					OddLineSegU2.push_back(se);
				}
			}
		}

		std::vector<SimEdge> OddLineSegE1;
		if(OddPointsE1.size() > 0)
		{
			for(i=0; i<OddPointsE1.size() - 1; i++)
			{
				for(j=i+1; j<OddPointsE1.size(); j++)
				{
					SimEdge se;	 se.A = OddPointsE1[i].Ver; se.B = OddPointsE1[j].Ver; se.x = 0; se.y = 0; se.z = 0;
					OddLineSegE1.push_back(se);
				}
			}
		}

		/* ����Ƿ��б��������ɵı����档*/
		if((OddLineSegU2.size() == 0) && (OddLineSegE1.size() == 0)) break;

		if(OddLineSegU2.size() > 0)
		{
			std::vector<SimEdge> remain;
			for( std::vector<SimEdge>::iterator p = kept.begin(); p != kept.end(); p++ )
			{
				int iFlag = 1;
				for( std::vector<SimEdge>::iterator q = OddLineSegU2.begin(); q != OddLineSegU2.end(); q++ )
				{
					if( *p == *q )	{ hold.push_back(*p);	iFlag = 0; break;	}
				}
				if(iFlag)	remain.push_back(*p);
			}
	
			kept.clear();
			for( std::vector<SimEdge>::iterator p = remain.begin(); p != remain.end(); p++ )
				kept.push_back(*p);
		}
		
		if(OddLineSegE1.size() > 0)
		{
			for( std::vector<SimEdge>::iterator p = hold.begin(); p != hold.end(); p++ )
			{
				for( std::vector<SimEdge>::iterator q = OddLineSegE1.begin(); q != OddLineSegE1.end(); q++ )
				{
					if( *p == *q )	{ kept.push_back(*p); break;	}
				}
			}
		}
		
		if(lastEdgeSize == kept.size()) break;
		lastEdgeSize = kept.size();
	}

#if 0
	sprintf(debugFile, "%s\\Flt2.ive", tempDirectory);
	OutPutEdgesIve(kept, debugFile);
#endif

	/* ����LINE_LOOP�������������ı�����һ����LINE_LOOP*/
	std::vector<SimEdge> alledges;	alledges.clear();
	for( std::vector<SimEdge>::iterator q = kept.begin(); q != kept.end(); q++ )  alledges.push_back(*q);
	osg::ref_ptr<osg::Vec3dArray> perimeter = new osg::Vec3dArray;
	while(alledges.size() > 0)
	{
    	/* �����һ��LINE_LOOP */
		osg::ref_ptr<osg::Vec3dArray> coords = new osg::Vec3dArray;
		////////////////////////////////////////////////////////////////////
		SimEdge e = alledges.front();
		osg::Vec3d A = e.A;
		osg::Vec3d B = e.B;
		unsigned int iNum = 0;
		int iFlag = 1;
		do {
			for( std::vector<SimEdge>::iterator p = alledges.begin(); p != alledges.end(); p++ )
			{
				if( p->has(B) && !(e == *p) )
				{
					coords->push_back( A );
					A = B;
					if( B == p->A )
						B = p->B;
					else
						B = p->A;
	
					e = *p;
					break;
				}
			}
			iNum++;
			if(iNum > alledges.size()) {iFlag = 0; break; }
		}while( !(e == alledges.front()) );
		coords->push_back(A);

		/* ��perimeter��LINE_LOOP���뵽һ��geometry���� */
		if (coords->size()>=3)
		{
			int iFlag = 0;
			/* �ж�perimeter��������㣬���а����м�ֵ�㣬˵�����LineLoop���Ǳ�Ե��LineLoop��*/
			for( osg::Vec3dArray::iterator pvx = coords->begin(); pvx != coords->end(); pvx++)
			{
				osg::Vec3d v; v = *pvx;
				if(v.x() == mNode.maxX) iFlag++;
				if(v.y() == mNode.maxY) iFlag++;
				if(v.x() == mNode.minX) iFlag++;
				if(v.y() == mNode.minY) iFlag++;
			}
			
			if(iFlag > 0)
			{
				for( osg::Vec3dArray::iterator pvx = coords->begin(); pvx != coords->end(); pvx++) 	perimeter->push_back(*pvx);
				//break;	
			}
		}

		/* ������keep����ԭΪtrue */
		for( std::vector<SimEdge>::iterator p = alledges.begin(); p != alledges.end(); p++ ) 
		p->keep = true;
		
		/* ��alledges�����߼��������perimeters */
		for( osg::Vec3dArray::iterator v = coords->begin(); v != coords->end(); v++ )
		{
			for( std::vector<SimEdge>::iterator e = alledges.begin(); e != alledges.end(); e++ )
			{
				if(e->has(*v)) e->keep = false;
			}
		}

		/* ��keep == false�Ķ�ȥ�� */
		std::vector<SimEdge> kept;
		for( std::vector<SimEdge>::iterator p = alledges.begin(); p != alledges.end(); p++ )
		{
			if( p->keep )
				kept.push_back( *p );
		}
		alledges.clear();				    
		for( std::vector<SimEdge>::iterator e = kept.begin(); e != kept.end(); e++ )
		{
			alledges.push_back( *e );
		}
		
	}	//while

#if 0
	sprintf(debugFile, "%s\\Flt3.ive", tempDirectory);
	OutPutLineLoopIve2(perimeter, debugFile);
#endif

	/* ���ҵ���Line_Loop�����һ��Geometry�� */
	/* �������������е�����һ����LINE_LOOP*/
	if(perimeter->size() == 0) return 0;

	osg::ref_ptr<osg::Vec3Array> vx = new osg::Vec3Array;
	for( osg::Vec3dArray::iterator pvx = perimeter->begin(); pvx != perimeter->end(); pvx++) 
	{
		osg::Vec3 v = *pvx;
		vx->push_back(v);
	}
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	gnode->addDrawable(gm);
	gm->setVertexArray( vx.get() );
	gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, vx->size()) );

#endif	//USE_TOTALEDGES
	if(!(gnode.valid())) return 0;
	return 1;
}


/*******************************************************
 *
 *  ���IVE�ļ�����ز����� (�����ô��룬�ǹ���ʹ�á�)
 *
 *******************************************************/

/* ��ȡһ��IVE�ļ��������ж������Ϊ��γ�߸�ʽ�����Բ鿴�м����õĺ��� */
void ReadIVEAndChangeVertexToGeodesy(char *fname)
{
	FILE *fp;
	int idx = 0;
	Geodesy geod;
	Carton  cart;

	fp = fopen("E:\\zbaavx.txt", "wt"); 
	osg::ref_ptr<osg::Node> ive = osgDB::readNodeFile(fname);
	osg::ref_ptr<osg::Group> gp	= dynamic_cast<osg::Group*> (ive.get());
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));		
	    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    	{
	        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
	        if( geometry.valid() )
	        {
				osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
				for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
				{
					osg::Vec3 v; v = *pvx;
					cart.x = v[0]; cart.y = v[1]; cart.z = v[2];
					CartToGeod(&cart, &geod);
					fprintf(fp, "%3d: %12.6f  %12.6f  %12.6f\n", idx++, geod._lon, geod._lat, geod._elevation);
				}
			}
		}
	}
	fclose(fp);
}

/* ��ȡһ��IVE�ļ������������һ�����LINE_LOOP�����Բ鿴�м����õĺ��� */
void ReadIVEAndComputePerimeter(char *fname)
{
	int idx = 0;

	osg::ref_ptr<osg::Node> ive = osgDB::readNodeFile(fname);
	osg::ref_ptr<osg::Group> gp	= dynamic_cast<osg::Group*> (ive.get());
    std::vector<SimEdge>  totalEdges;	totalEdges.clear();
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
		
		if(("Grass" == gnode->getName().substr(0, strlen("Grass"))) || ("pc_yellow" == gnode->getName().substr(0, strlen("pc_yellow"))))
			continue;
		
	    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    	{
	        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
			osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray() );
	        if( geometry.valid() )
	        {
			    for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
			    {
			        osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
			        
					switch( pset->getType() )
					{
					    case osg::PrimitiveSet::DrawArraysPrimitiveType:
					    {
					        osg::DrawArrays *da = dynamic_cast<osg::DrawArrays *>(pset.get());
					        int first = da->getFirst();
					        unsigned int count = da->getCount();
					
					        switch(pset->getMode() )
					        {
					            case osg::PrimitiveSet::TRIANGLES:
					            {
									for( unsigned int ii = 0; ii < count; ii+=3 )
									{
										osg::Vec3 vx, vy, vz;
										vx = (*coords)[first + ii + 0];
										vy = (*coords)[first + ii + 1];
										vz = (*coords)[first + ii + 2];
										SimEdge se;	 
										if( vx != vy)
										{
											se.A = vx; se.B = vy; se.x = 0; se.y = 0; se.z = 0;
											totalEdges.push_back(se);
										}
										if( vy != vz)
										{
											se.A = vy; se.B = vz; se.x = 0; se.y = 0; se.z = 0;
											totalEdges.push_back(se);
										}
										if( vx != vz)
										{
											se.A = vx; se.B = vz; se.x = 0; se.y = 0; se.z = 0;
											totalEdges.push_back(se);
										}
									}
					            }
					            break;
					        }
					    }
					    break;
					}
			    }
			}
		}
	}
	
#if 1
	/* �ڶ������ж�ÿһ����Ч�ߣ��Ƿ��Ǳ߽�ıߣ������Ǵ�ֱƽ�����ϵ����������㣬������ͽ����ཻ���Ͳ��Ǳ߽�� */
	std::vector<SimEdge> peri; peri.clear();
	int i = 0;
    for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
    {
    	printf("%d\n", i++);

		SimEdge se = *p;
    	if(DetectIntersectOfEdge(&se, ive) == false)
    		peri.push_back(se);
    }

	/* �����е�Edge�����ֻ�е���ͬ���Edge���ӣ�����������Edge*/
	totalEdges.clear();
    for( std::vector<SimEdge>::iterator p = peri.begin(); p != peri.end(); p++ )
	{
		SimEdge se = *p;
		totalEdges.push_back(se);
	}

	OutPutEdgesIve(peri, "E:\\Zbaa_LineLoop.ive");

#else
	/* �����һ��LINE_LOOP */
	osg::ref_ptr<osg::Vec3Array> perimeter = new osg::Vec3Array;
	if(computePerimeter( totalEdges, *perimeter.get()) == 0) return;
		
	OutPutLineLoopIve(perimeter, "E:\\Zbaa_LineLoop.ive");
#endif

}

/* ��ȡһ��IVE�ļ���ɾ�����˵��ߣ����Բ鿴�м����õĺ��� */
void ReadIVEAndDeleteLine(char *fname)
{
	int idx = 0;

	osg::ref_ptr<osg::Node> ive = osgDB::readNodeFile(fname);
	osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (ive.get());
    std::vector<SimEdge>  totalEdges;	totalEdges.clear();
    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
	{
        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
		osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray() );
        if( geometry.valid() )
        {
		    for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
		    {
		        osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
		        
				switch( pset->getType() )
				{
				    case osg::PrimitiveSet::DrawArraysPrimitiveType:
				    {
				        osg::DrawArrays *da = dynamic_cast<osg::DrawArrays *>(pset.get());
				        int first = da->getFirst();
				        unsigned int count = da->getCount();
				
				        switch(pset->getMode() )
				        {
				            case osg::PrimitiveSet::LINES:
				            {
								for( unsigned int ii = 0; ii < count; ii+=2 )
								{
									osg::Vec3 vx, vy;
									vx = (*coords)[first + ii + 0];
									vy = (*coords)[first + ii + 1];
									SimEdge se;	 
									if( vx != vy)
									{
										se.A = vx; se.B = vy; se.x = 0; se.y = 0; se.z = 0;
										totalEdges.push_back(se);
									}
								}
				            }
				            break;
				        }
				    }
				    break;
				}
		    }
		}
	}

	/* �����е�Edge�����ֻ�е���ͬ���Edge���ӣ�����������Edge*/
	std::vector<SimEdge> peri; peri.clear();
    for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
	{
		int iFlagA, iFlagB;
		iFlagA = iFlagB = 0;

		osg::Vec3 A = p->A;
		osg::Vec3 B = p->B;

	    for( std::vector<SimEdge>::iterator q = totalEdges.begin(); q != totalEdges.end(); q++ )
	    {
	    	if ((!( *p == *q )) && ( q->has(A)))	iFlagA = 1;
	    	if ((!( *p == *q )) && ( q->has(B)))	iFlagB = 1;
	    }
		if((iFlagA) && (iFlagB))
		{
			SimEdge se = *p;
			peri.push_back(se);
		}
	}

	OutPutEdgesIve(peri, "E:\\Zbaa_LineLoop2.ive");
}

/* ��ȡһ��IVE�ļ������������һ�����LINE_LOOP�����ive�ļ���fltת�������Բ鿴�м����õĺ���*/
void ReadIVEAndConnectLineLoop(char *fname)
{
	int idx = 0;

	osg::ref_ptr<osg::Node> ive = osgDB::readNodeFile(fname);
	osg::ref_ptr<osg::Group> gp	= dynamic_cast<osg::Group*> (ive.get());
    std::vector<SimEdge>  totalEdges;	totalEdges.clear();
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
		
	    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    	{
	        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
			osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray() );
	        if( geometry.valid() )
	        {
			    for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
			    {
			        osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
			        
					switch( pset->getType() )
					{
					    case osg::PrimitiveSet::DrawArraysPrimitiveType:
					    {
					        osg::DrawArrays *da = dynamic_cast<osg::DrawArrays *>(pset.get());
					        int first = da->getFirst();
					        unsigned int count = da->getCount();
					        switch(pset->getMode() )
					        {
					            case osg::PrimitiveSet::LINE_STRIP:
					            {
									for( unsigned int ii = 1; ii < count; ii++ )
									{
										osg::Vec3 vx, vy;
										vx = (*coords)[first + ii - 1];
										vy = (*coords)[first + ii - 0];
										SimEdge se;	 
										if( vx != vy)
										{
											se.A = vx; se.B = vy; se.x = 0; se.y = 0; se.z = 0; se.keep = true;
											totalEdges.push_back(se);
										}
									}
					            }
					            break;
					            case osg::PrimitiveSet::LINES:
					            {
									for( unsigned int ii = 0; ii < count; ii+=2 )
									{
										osg::Vec3 vx, vy;
										vx = (*coords)[first + ii + 0];
										vy = (*coords)[first + ii + 1];
										SimEdge se;	 
										if( vx != vy)
										{
											se.A = vx; se.B = vy; se.x = 0; se.y = 0; se.z = 0; se.keep = true;
											totalEdges.push_back(se);
										}
									}
					            }
					            break;
					        }
					    }
					    break;
					}
			    }
			}
		}
	}
	
	/* ���ҳ���ͷ�����ж��ٸ���ͷ�ĵ� */
	std::vector<osg::Vec3> VxNc; VxNc.clear();
    for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
	{
		int iFlagA, iFlagB;
		iFlagA = iFlagB = 0;

		osg::Vec3 A = p->A;
		osg::Vec3 B = p->B;

	    for( std::vector<SimEdge>::iterator q = totalEdges.begin(); q != totalEdges.end(); q++ )
	    {
	    	if ((!( *p == *q )) && ( q->has(A)))	iFlagA = 1;
	    	if ((!( *p == *q )) && ( q->has(B)))	iFlagB = 1;
	    }
		if(!iFlagA)
		{
			VxNc.push_back(A);
			p->keep = false;
		}
		if(!iFlagB)
		{
			VxNc.push_back(B);
			p->keep = false;
		}
	}

	std::vector<SimEdge> peri; peri.clear();
    for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
	{
		if(p->keep)
			peri.push_back(*p);
	}

//	OutPutEdgesIve(peri, "E:\\gm.ive");
	totalEdges.clear();
    for( std::vector<SimEdge>::iterator p = peri.begin(); p != peri.end(); p++ )
		totalEdges.push_back(*p);

	/* �����е�����һ����LINE_LOOP*/
	osg::ref_ptr<osg::Vec3Array> coords = new osg::Vec3Array;
	SimEdge e = totalEdges.front();
	osg::Vec3 A = e.A;
	osg::Vec3 B = e.B;
	do {
		for( std::vector<SimEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
		{
			if( p->has(B) && !(e == *p) )
			{
				coords->push_back( A );
				A = B;
				if( B == p->A )
					B = p->B;
				else
					B = p->A;

				e = *p;
				break;
			}
		}
	}while( !(e == totalEdges.front()) );
	coords->push_back(A);

	OutPutLineLoopIve(coords, "E:\\genz.ive");

}

/*******************************************************
 *
 *  ����һ��������ţ��ҵ�����������ڵ��� Group��btg������ (�����ô��룬�ǹ���ʹ�á�)
 *
 *******************************************************/

typedef struct _VertValid_
{
	int index;
	bool isOK;
}VertValid;

typedef struct _OneTriangle_
{
	int A, B, C;
}OneTriangle;

std::vector<VertValid>  findVertice;  	/* ��Ч�㼯�ϣ��������Ŷ���ֱ�ӻ��߼�����������ж���ļ��� */
std::vector<OneTriangle> tris;			/* ��Ч�����Σ�������Ч������������μ���������ҳ����� */
	
/* �ж�һ����������Ƿ�����Ч�㼯�����棬���Բ鿴�м����õĺ��� */
bool CheckVertexIndexInFV(int idx)
{
	bool ret = false;
	for( std::vector<VertValid>::iterator m = findVertice.begin(); m != findVertice.end(); m++ )
	{
		if(m->index == idx) { ret = true; break;}
	}
	return ret;	
}

/* ����һ��������ţ��ҵ�����������ڵ��� Group��btg�����������Բ鿴�м����õĺ��� */
void FindFaceWithSpecVertex(char *btgName, int idx)
{
	/* ���ȶ���btg�ļ����������һ��Group */
	ReadBTGToDataStruct(btgName);
		
	/* �ȼ�������еı�index����, ������totalEdges���� */
    std::vector<IndEdge>  totalEdges;	totalEdges.clear();
	for(unsigned int i=0;i<curTri.size();i++)
	{
		for(unsigned int j=0;j<curTri[i].mElement.size();j++)
		{
			for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size(); k+=3)
			{
				int a, b, c;
				a = curTri[i].mElement[j].mVertex[k + 0].P;
				b = curTri[i].mElement[j].mVertex[k + 1].P;
				c = curTri[i].mElement[j].mVertex[k + 2].P;
				IndEdge se;
				se.A = a; se.B = b; se.keep = true; totalEdges.push_back(se);
				se.A = b; se.B = c; se.keep = true; totalEdges.push_back(se);
				se.A = c; se.B = a; se.keep = true; totalEdges.push_back(se);
			}//k
		}// j
	}// i

	/* �ҳ����ض����㹲�ߵ����е㣬������ findVertice ���� */
	findVertice.clear();
	VertValid first;  first.index = idx; first.isOK = false; findVertice.push_back(first);

	int lastIdx = -1;
	while(1)
	{	
		int curIdx, exitFlag;

		exitFlag = 1;
		for( std::vector<VertValid>::iterator q = findVertice.begin(); q != findVertice.end(); q++ )
		{
			if(q->isOK == false) 
			{ 
				curIdx = q->index;
				if(lastIdx == curIdx)
				{
					q->isOK = true;
					continue;
				}
				exitFlag = 0; 
				break; 
			}
		}
		if(exitFlag) 
			break;
		
		for( std::vector<IndEdge>::iterator p = totalEdges.begin(); p != totalEdges.end(); p++ )
		{
			if(p->has(curIdx))
			{
				int next, iFlag = 0;
				if(curIdx == p->A) next = p->B; else next = p->A;
				for( std::vector<VertValid>::iterator m = findVertice.begin(); m != findVertice.end(); m++ )
				{
					if(m->index == curIdx) {m->isOK = true;}
					if(m->index == next) {iFlag = 1; break;}
				}
				if(!iFlag)
				{
					VertValid vv; vv.index = next; vv.isOK = false; findVertice.push_back(vv);
				}
			}
		}//p
		
		lastIdx = curIdx;
	}//while
	
	/* ������Ч�㼯�� findVertice �������е���Ч�����Σ�ֻҪ�����ε�һ����������Ч�㣬�ͽ���������μ�¼�� tris ���� */
	tris.clear();
	for(unsigned int i=0;i<curTri.size();i++)
	{
		for(unsigned int j=0;j<curTri[i].mElement.size();j++)
		{
			for(unsigned int k=0;k<curTri[i].mElement[j].mVertex.size(); k+=3)
			{
				int a, b, c;
				a = curTri[i].mElement[j].mVertex[k + 0].P;
				b = curTri[i].mElement[j].mVertex[k + 1].P;
				c = curTri[i].mElement[j].mVertex[k + 2].P;
				
				/* ֻ���ж�һ���㣬������������Ȼ�Ͷ��������ˡ�*/
				if(CheckVertexIndexInFV(a))
				{
					OneTriangle ot; ot.A = a; ot.B = b; ot.C = c; tris.push_back(ot);
				}
			}//k
		}// j
	}// i	
}

/* ��һ����ָ�������μ��뵽һ��ָ�����ʵ������μ������棬���Բ鿴�м����õĺ��� */
int AddTrianglesToCurTri(char *mater)
{
	for(unsigned int i=0;i<curTri.size();i++)
	{
		curTri[i].mElement.clear();
		if(!strcmp(mater, curTri[i].cMtype))
		{
			TriElement te; te.mVertex.clear();
			for(unsigned int j=0; j<tris.size(); j++)
			{
				TriVertex tv;
				tv.P = tris[j].A; tv.N = 0; tv.T = tris[j].A; te.mVertex.push_back(tv);
				tv.P = tris[j].B; tv.N = 0; tv.T = tris[j].B; te.mVertex.push_back(tv);
				tv.P = tris[j].C; tv.N = 0; tv.T = tris[j].C; te.mVertex.push_back(tv);                          
			}                                                                         
			te.iSum = te.mVertex.size();
			curTri[i].mElement.push_back(te);
		}
	}
	
	return 1;
}

/* ���һ�������ε������㣬�Ƿ���һ������tris�����μ������棬���Բ鿴�м����õĺ��� */
int CheckTri(int a, int b, int c)
{
	int ret = 0;
	for(unsigned int i=0; i<tris.size(); i++)
	{
		if( ((a==tris[i].A) || (a==tris[i].B) || (a==tris[i].C)) &&
			((b==tris[i].A) || (b==tris[i].B) || (b==tris[i].C)) &&
			((c==tris[i].A) || (c==tris[i].B) || (c==tris[i].C)) )
		{ ret = 1; break; }
	}
	return ret;
}

/* ��һ����ָ�������δ�һ��ָ�����ʵ������μ��������Ƴ������Բ鿴�м����õĺ��� */
int RemoveTrianglesFromCurTri(char *mater)
{
	int indexType;
	unsigned int k;
	
	for(unsigned int i=0;i<curTri.size();i++)
	{
		if(!strcmp(mater, curTri[i].cMtype))
		{
			indexType = curTri[i].cPropVal;

			for(unsigned int j=0;j<curTri[i].mElement.size();j++)
			{
				/* �ҵ�ƥ�������β���ֵΪ -1 */
				for(k=0;k<curTri[i].mElement[j].mVertex.size();k+=3)
				{
					int a, b, c;
					a = curTri[i].mElement[j].mVertex[k + 0].P;
					b = curTri[i].mElement[j].mVertex[k + 1].P;
					c = curTri[i].mElement[j].mVertex[k + 2].P;
					if(CheckTri(a, b, c))
					{
						curTri[i].mElement[j].mVertex[k + 0].P = -1;
						curTri[i].mElement[j].mVertex[k + 1].P = -1;
						curTri[i].mElement[j].mVertex[k + 2].P = -1;
					}
				}
				
				/* ��������ǰ*/
				k = 0;
				std::vector<TriVertex> tvlist; tvlist.clear();
				while(k < curTri[i].mElement[j].mVertex.size())
				{
					if(curTri[i].mElement[j].mVertex[k].P != -1)
					{
						TriVertex tv;
						tv.P = curTri[i].mElement[j].mVertex[k].P; tv.N = curTri[i].mElement[j].mVertex[k].N; tv.T = curTri[i].mElement[j].mVertex[k].T; tvlist.push_back(tv);
					}
					k++;
				}
				curTri[i].mElement[j].mVertex.clear();
				for(std::vector<TriVertex>::iterator p = tvlist.begin(); p != tvlist.end(); p++)
					curTri[i].mElement[j].mVertex.push_back(*p);
				curTri[i].mElement[j].iSum = curTri[i].mElement[j].mVertex.size();
			}
		}
	}
	
	return 1;
}

/* �޸�btg�ļ������������Index��Ŷ������һ����Ĳ��ʴ�Stream����ת��ΪEvergreenForest���ͣ����Բ鿴�м����õĺ��� */
void ModifyBtg(char *btgname, int Index)
{
	FindFaceWithSpecVertex(btgname, Index);
	AddTrianglesToCurTri("Stream");
	RemoveTrianglesFromCurTri("EvergreenForest");
}




/*******************************************************
 *
 *  ����һ�����ǻ��˵�Group�� �����ڲ���������㣬����������ƽ���ϵ�ɭ�֡�
 *
 *******************************************************/
#define WOOD_COVERAGE	10000.0
#define TREE_DENSITY  1.0
#define WOOD_SIZE 0.0

int idxOfRandomTree;
double length(const osg::Vec3 v)
{ return sqrt(v.x() * v.x() + v.y() * v.y() + v.z() * v.z()); }

osg::Vec3 cross(const osg::Vec3 v1, const osg::Vec3 v2)
{
  return osg::Vec3(v1.y()*v2.z() - v1.z()*v2.y(),
                   v1.z()*v2.x() - v1.x()*v2.z(),
                   v1.x()*v2.y() - v1.y()*v2.x());
}
osg::Vec3 normalize(const osg::Vec3 v)
{
  double len = 1.0 / length(v);
  osg::Vec3 norm; norm.set( len * v.x(), len * v.y(), len * v.z());
  return norm;
}

/* ����õ��������� */
string getRandomTreeName()
{
	int idx = int( mt_rand(&random_seed) * 10.0);
	char treeName[12]; sprintf(treeName, "Forest_%1d", idx);
	return treeName;
}

string getTreeName(string curr)
{
	/* �ȼ�鵱ǰ�����Ƿ�Ϸ� */
	string currName = curr;
	splitMaterialName(currName);
	if((materialParts[0] == "Forest") || (materialParts[0] == "forest") || (materialParts[0] == "Tree") || (materialParts[0] == "tree"))
	{
		int id;
		if(materialParts.size() > 1)
		{
			id = atoi(materialParts[1].c_str());
			if((id>=0) && (id <= 99))
			{
				char treeName[12]; sprintf(treeName, "Forest_%1d", id);
				currName = treeName;
			}else
				currName = getRandomTreeName();
		}else
			currName = getRandomTreeName();
	}else
		currName = getRandomTreeName();
	return currName;
}

/* ��֪������ĳһֱ�������v������õ㾭γ�ȣ���������õ�߶ȣ���ת��Ϊֱ������㲢���ء�*/
void getActualVertex(osg::Vec3 v, int mode, string grName)
{
	Geodesy geod;
	Carton  cart;
	double  actualElev;
	
	cart.x = v.x(); cart.y = v.y(); cart.z = v.z();
	CartToGeod(&cart, &geod);
	geod._elevation = actualElev = GetElevation(geod._lon, geod._lat, 1.0);
	GeodToCart(&geod, &cart);

	char tile[16];
	Bucket bt;
	set_bucket(&bt, geod._lon, geod._lat);
	gen_index_str(&bt, tile);
	
	SimPoint sp; sp.x  = cart.x; sp.y  = cart.y; sp.z  = cart.z; sp.lon = geod._lon; sp.lat = geod._lat; sp.elev = geod._elevation; sp.index = idxOfRandomTree++; sp.tile = tile; 
	if(mode == 0)
	{
		sp.name = getTreeName(grName);
#if 1
		/* �������򣬰�������һ��߶ȡ� */
		geod._elevation -= 1.5;
		GeodToCart(&geod, &cart);
		sp.x  = cart.x; sp.y  = cart.y; sp.z  = cart.z; sp.elev = geod._elevation;
#endif
		TreePoints.push_back(sp);
	}
	if(mode == 1)
	{
		sp.name = "USR_WHITE_LIGHTS_0_30000_" + sp.tile;
		CultPoints.push_back(sp);
	}
}

/* ��һ����������������������㣬��Ϊ�������������ã������������ǵƣ���mode������ mode = 0--����mode = 1--�� */
void CreateRandomPointsInTrangle(osg::Vec3 v0, osg::Vec3 v1, osg::Vec3 v2, int mode, string grName) 
{
	float wood_size = WOOD_SIZE;
	float tree_density = TREE_DENSITY;
	float wood_coverage = mWoodCoverage;
	osg::Vec3 normal = cross(v1 - v0, v2 - v0);

	// Compute the area
	float area = 0.5f*length(normal);
	if (area <= 10.0)
		return;

	// For partial units of area, use a zombie door method to
	// create the proper random chance of a point being created
	// for this triangle
	float unit = area + mt_rand(&random_seed)*wood_coverage;
	int woodcount = (int) (unit / wood_coverage);
	
	for (int j = 0; j < woodcount; j++) {
	
		if (wood_size < area) {
			// We need to place a wood within the triangle and populate it
			
			// Determine the center of the wood
			float x = mt_rand(&random_seed);
			float y = mt_rand(&random_seed);
			
			// Determine the size of this wood in m^2, and the number
			// of trees in the wood
			float ws = wood_size + wood_size * (mt_rand(&random_seed) - 0.5f);
			unsigned total_trees = ws / tree_density;
			float wood_length = sqrt(ws);
			
			// From our wood size, work out the fraction on the two axis.
			// This will be used as a factor when placing trees in the wood.
			float x_tree_factor = wood_length / length(v1 -v0);
			float y_tree_factor = wood_length / length(v2 -v0);
			
			for (unsigned k = 0; k <= total_trees; k++) {
			
				float a = x + x_tree_factor * (mt_rand(&random_seed) - 0.5f);
				float b = y + y_tree_factor * (mt_rand(&random_seed) - 0.5f);
				
				
				// In some cases, the triangle side lengths are so small that the
				// tree_factors become so large as to make placing the tree within
				// the triangle almost impossible. In this case, we place them
				// randomly across the triangle.
				if (a < 0.0f || a > 1.0f) a = mt_rand(&random_seed);
				if (b < 0.0f || b > 1.0f) b = mt_rand(&random_seed);
				
				if ( a + b > 1.0f ) {
					a = 1.0f - a;
					b = 1.0f - b;
				}
				
				float c = 1.0f - a - b;
				
				osg::Vec3 randomPoint;
				randomPoint.x() = a*v0.x() + b*v1.x() + c*v2.x();
				randomPoint.y() = a*v0.y() + b*v1.y() + c*v2.y();
				randomPoint.z() = a*v0.z() + b*v1.z() + c*v2.z();
				
				getActualVertex(randomPoint, mode, grName);
			}
		} else {
			// This triangle is too small to contain a complete wood, so just
			// distribute trees across it.
			unsigned total_trees = area / tree_density;
			
			for (unsigned k = 0; k <= total_trees; k++) {
			
				float a = mt_rand(&random_seed);
				float b = mt_rand(&random_seed);
				
				if ( a + b > 1.0f ) {
					a = 1.0f - a;
					b = 1.0f - b;
				}
				
				float c = 1.0f - a - b;
				
				osg::Vec3 randomPoint;
				randomPoint.x() = a*v0.x() + b*v1.x() + c*v2.x();
				randomPoint.y() = a*v0.y() + b*v1.y() + c*v2.y();
				randomPoint.z() = a*v0.z() + b*v1.z() + c*v2.z();
				
				getActualVertex(randomPoint, mode, grName);
			}
		}
	}
}
 
/* ��һ���������������ɴ��������㣬��Ϊ���ɭ���á�*/
void CreateTreePointsFromGroup(osg::ref_ptr<osg::Group> frst)
{
	CreateRandomPointsFromGroup(frst, 0);
}

/* ��һ���������������ɴ��������㣬��Ϊ������á�*/
void CreateCultPointsFromGroup(osg::ref_ptr<osg::Group> lights)
{
	CreateRandomPointsFromGroup(lights, 1);
}

/* ��һ���������������ɴ��������㣬��Ϊ�������������ã������������ǵƣ���mode������ mode = 0--����mode = 1--�� */
void CreateRandomPointsFromGroup(osg::ref_ptr<osg::Group> frst, int mode)
{
	unsigned int i;
	string grName;
	if(mode == 0) TreePoints.clear();
	if(mode == 1) CultPoints.clear();
	idxOfRandomTree = 0;
	
	/* �ҳ���������㣬*/
	for(i = 0; i<frst->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (frst->getChild(i));
		grName = gnode->getName();
		for(unsigned int j = 0; j<gnode->getNumDrawables(); j++)
		{
			osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(j));
			osg::Vec3Array* coords = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
	
			for (unsigned int ipr=0; ipr< geometry->getNumPrimitiveSets(); ipr++) 
			{
				osg::ref_ptr<osg::PrimitiveSet> prset = geometry->getPrimitiveSet(ipr);
				switch( prset->getType() )
				{
					case osg::PrimitiveSet::DrawArraysPrimitiveType:
						{
							osg::DrawArrays* gdui = dynamic_cast<osg::DrawArrays*>(prset.get());
							if( gdui )
							{
								switch(prset->getMode() )
								{
									case osg::PrimitiveSet::TRIANGLES:
										for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
										{
											unsigned int a,b,c;
											a = gdui->index(i);
											b = gdui->index(i+1);
											c = gdui->index(i+2);	
											
											CreateRandomPointsInTrangle(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), mode, grName);
										}
										break;
									default:
										break;
								}
							}
						 }
						break;
					case osg::PrimitiveSet::DrawArrayLengthsPrimitiveType:
					case osg::PrimitiveSet::DrawElementsUBytePrimitiveType:
					case osg::PrimitiveSet::DrawElementsUShortPrimitiveType:
						break;
					case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
						{
							osg::DrawElementsUInt *gdui = dynamic_cast<osg::DrawElementsUInt *>(prset.get());
							if( gdui )
							{
								switch(prset->getMode() )
								{
									case osg::PrimitiveSet::TRIANGLES:
										for(unsigned int i=0; i<gdui->getNumIndices(); i+=3)
										{
											unsigned int a,b,c;
											a = gdui->index(i);
											b = gdui->index(i+1);
											c = gdui->index(i+2);	

											CreateRandomPointsInTrangle(*(coords->begin()+a), *(coords->begin()+b), *(coords->begin()+c), mode, grName);
										}
										break;
									default:
										break;
								}
							}
						}
						break;
					default:
						break;
				}	//case
			}	//for ipr
		}	//for j
	} // for i
}

/* ��TreePoints���������name���컯��ʵ�ֲ�ͬ����ֱ���ʾ��Ч����*/
int AddLODIntoTreePoints()
{
	if(TreePoints.size() > 0)
	{
		for(unsigned int i=0; i<TreePoints.size(); i++)
		{
			SimPoint *sp = &TreePoints[i];
			
			/* �ȼ�鵱ǰ�����Ƿ�Ϸ� */
			string currName = getTreeName(sp->name);
			
			if((i%8) == 0)	sp->name = currName + "_0_0_10000_" + sp->tile;
			if((i%8) == 1)	sp->name = currName + "_0_0_10000_" + sp->tile;
			if((i%8) == 2)	sp->name = currName + "_0_0_6000_" + sp->tile;
			if((i%8) == 3)	sp->name = currName + "_0_0_6000_" + sp->tile;
			if((i%8) == 4)	sp->name = currName + "_0_0_3000_" + sp->tile;
			if((i%8) == 5)	sp->name = currName + "_0_0_3000_" + sp->tile;
			if((i%8) == 6)	sp->name = currName + "_0_0_3000_" + sp->tile;
			if((i%8) == 7)	sp->name = currName + "_0_0_3000_" + sp->tile;
		}
	}
	return 1;
}

/* ��TreePoints�������������ת��ΪBtg��ʽ�����ݡ�*/
int CreateTreesDataStructFromTreePoints()
{	
	char sMater[64];

	//��ʼ���ڴ�
	curPnt.clear();		totalNormal = 0;
	curText.clear();
	curTri.clear();

	mNode.maxX = mNode.maxY = mNode.maxZ = -FLT_MAX;	mNode.minX = mNode.minY = mNode.minZ = FLT_MAX;
	if(TreePoints.size() > 0)
	{
		unsigned int i;
		std::vector<string> mats;	mats.clear();
		
		/* д������Ϣ */
		for(i=0; i<TreePoints.size(); i++)
		{
			TerPoint tp; tp.px = TreePoints[i].x; tp.py = TreePoints[i].y; tp.pz = TreePoints[i].z; tp.nx = 0; tp.ny = 0; tp.nz = 128; curPnt.push_back(tp);
			TreePoints[i].index = i;
			
			/*��ֵ������mNode ����*/
			if(mNode.maxX < TreePoints[i].x) mNode.maxX = TreePoints[i].x;
			if(mNode.maxY < TreePoints[i].y) mNode.maxY = TreePoints[i].y;
			if(mNode.maxZ < TreePoints[i].z) mNode.maxZ = TreePoints[i].z;
			if(mNode.minX > TreePoints[i].x) mNode.minX = TreePoints[i].x;
			if(mNode.minY > TreePoints[i].y) mNode.minY = TreePoints[i].y;
			if(mNode.minZ > TreePoints[i].z) mNode.minZ = TreePoints[i].z;

			/* ͳ��һ���м��ֲ���*/
			int iFlag = 0;
			for(std::vector<string>::iterator p = mats.begin(); p != mats.end(); p++)
			{
				if(TreePoints[i].name == *p) { iFlag = 1; break;}
			}
			if(!iFlag) mats.push_back(TreePoints[i].name);
		}
	
		/* д�������� */
		totalNormal = 1;
		
		/* д����� */
		char cObjType;
		int indexType, indexOffset = 0;
		curTri.clear();
		cObjType = 9; indexType = 3;

		for(i = 0; i < mats.size(); i++)
		{
			TriType tt;
			tt.mElement.clear();

			/* ����te */
			TriElement te; te.mVertex.clear();
		    for( unsigned int n = 0; n < TreePoints.size(); n++ )
	    	{
	    		if(TreePoints[n].name == mats[i])
	    		{
					TriVertex tv;
					tv.P = tv.T = TreePoints[n].index; tv.N = 0; te.mVertex.push_back(tv);
	    		}
			}
			te.iSum = te.mVertex.size();
			tt.mElement.push_back(te);
			tt.iElems = tt.mElement.size();

			/* д���� */
			sprintf(sMater, "%s", mats[i].c_str());
			sprintf(tt.cMtype, "%s", sMater);
			
			tt.cObjType = cObjType;
			tt.cPropVal = indexType;
			curTri.push_back(tt);
		}
	
		TreePoints.clear();
	
		/* дbtg�ļ�ͷ */
		WriteBTGHeader();
	}
	return 1;
}

/* дBTG�ļ�ͷ��*/
void WriteBTGHeader()
{
	/* ����BTG�ļ�ͷ */
	iVersion = 6; iMagicNum = 0x5347; iCreationTime = 0x48170f2c; iNumOfToplvlObject = 5 + curTri.size();

	/* ��߽�������ĵ���Ϣ */
	dXpart = (mNode.maxX + mNode.minX) / 2.0;	
	dYpart = (mNode.maxY + mNode.minY) / 2.0;	
	dZpart = (mNode.maxZ + mNode.minZ) / 2.0;	
	fRadiusofBS = 0.0;
	Carton gbs_c, vert;
	double dist_squared, radius_squared = 0;
	gbs_c.x = dXpart; gbs_c.y = dYpart; gbs_c.z = dZpart;
	gbCenterX = dXpart; gbCenterY = dYpart; gbCenterZ = dZpart;
	
	/* �����Χ��İ뾶 */
	for(int j=0;j<(int)curPnt.size();j++)
	{
		vert.x = curPnt[j].px; vert.y = curPnt[j].py; vert.z = curPnt[j].pz;
		dist_squared = distance3Dsquared(&gbs_c, &vert);
		if ( dist_squared > radius_squared ) {
			radius_squared = dist_squared;
		}
	}
	fRadiusofBS = sqrt(radius_squared);
}

/*******************************************************
 *
 *  ����osg::Group��дһ��btg�ļ���
 *
 *******************************************************/

#define WRAPTIMES 100

class ArealConstraint: public  osgUtil::DelaunayConstraint 
{ 
public:
	int getinteriorTrisSize() { return _interiorTris.size(); }
};

/* ����һ��Group�ļ�ֵ */
void CalculateMaxMinOfGroup(osg::ref_ptr<osg::Group> gp)
{
	allFn.clear();
	FindMinMaxVisitor fmmv;
	gp->accept(fmmv);
	FindMaxMinInAllFn();
}

/* ��CultPoints���������name���컯��ʵ�ֲ�ͬ����ֱ���ʾ��Ч����*/
int AddLODIntoLightPoints()
{
	if(CultPoints.size() > 0)
	{
		for(unsigned int i=0; i<CultPoints.size(); i++)
		{
			SimPoint *sp = &CultPoints[i];
			GetWidthFromName(sp->name);
			string colorName = nameParts[0] + "_" + nameParts[1] + "_" + nameParts[2];
			if((i%4) == 0)	sp->name = colorName + "_0_50000_" + sp->tile;
			if((i%4) == 1)	sp->name = colorName + "_0_35000_" + sp->tile;
			if((i%4) == 2)	sp->name = colorName + "_0_20000_" + sp->tile;
			if((i%4) == 3)	sp->name = colorName + "_0_20000_" + sp->tile;
		}
	}
	return 1;
}


/* ����CultPoints�����ֵ������Ƶ�btg */
/* ��CultPoints�������������ת��ΪBtg��ʽ�����ݡ�*/
int CreateLightsDataStructFromCultPoints()
{
	char sMater[64];

	//��ʼ���ڴ�
	curPnt.clear();		totalNormal = 0;
	curText.clear();
	curTri.clear();

	mNode.maxX = mNode.maxY = mNode.maxZ = -FLT_MAX;	mNode.minX = mNode.minY = mNode.minZ = FLT_MAX;
	if(CultPoints.size() > 0)
	{
		unsigned int i;
		std::vector<string> mats;	mats.clear();
		
		/* д������Ϣ */
		for(i=0; i<CultPoints.size(); i++)
		{
			TerPoint tp; tp.px = CultPoints[i].x; tp.py = CultPoints[i].y; tp.pz = CultPoints[i].z; tp.nx = 0; tp.ny = 0; tp.nz = 128; curPnt.push_back(tp);
			CultPoints[i].index = i;
			
			/*��ֵ������mNode ����*/
			if(mNode.maxX < CultPoints[i].x) mNode.maxX = CultPoints[i].x;
			if(mNode.maxY < CultPoints[i].y) mNode.maxY = CultPoints[i].y;
			if(mNode.maxZ < CultPoints[i].z) mNode.maxZ = CultPoints[i].z;
			if(mNode.minX > CultPoints[i].x) mNode.minX = CultPoints[i].x;
			if(mNode.minY > CultPoints[i].y) mNode.minY = CultPoints[i].y;
			if(mNode.minZ > CultPoints[i].z) mNode.minZ = CultPoints[i].z;

			/* ͳ��һ���м��ֲ���*/
			int iFlag = 0;
			for(std::vector<string>::iterator p = mats.begin(); p != mats.end(); p++)
			{
				if(CultPoints[i].name == *p) { iFlag = 1; break;}
			}
			if(!iFlag) mats.push_back(CultPoints[i].name);
		}

		totalNormal = 1;

		/* д����� */
		char cObjType;
		int indexType, indexOffset = 0;
		curTri.clear();
		cObjType = 9; indexType = 3;

		for(i = 0; i < mats.size(); i++)
		{
			TriType tt;
			tt.mElement.clear();

			/* ����te */
			TriElement te; te.mVertex.clear();
		    for( unsigned int n = 0; n < CultPoints.size(); n++ )
	    	{
	    		if(CultPoints[n].name == mats[i])
	    		{
					TriVertex tv;
					tv.P = tv.T = CultPoints[n].index; tv.N = 0; te.mVertex.push_back(tv);
	    		}
			}
			te.iSum = te.mVertex.size();
			tt.mElement.push_back(te);
			tt.iElems = tt.mElement.size();

			/* д���� */
			sprintf(sMater, "%s", mats[i].c_str());
			sprintf(tt.cMtype, "%s", sMater);
			
			tt.cObjType = cObjType;
			tt.cPropVal = indexType;
			curTri.push_back(tt);
		}
		
		/* ����ת֮ǰ�����CultPints����ʾ�ƹ����ɭ���Ѿ�д�ꡣ*/
		CultPoints.clear();
		
		/* дbtg�ļ�ͷ */
		WriteBTGHeader();
	}
	return 1;
}

/* ��һ�������Ӷ��㼯�Ϻ��������꼯�� */
void CreateVecArrayAndAddTextureCoordinateIntoFace(osg::ref_ptr<osg::Group> gp)
{
	/* �������ֵ����Сֵ������ÿ��������������� */
	double Xzone = mNode.maxX - mNode.minX;
	double Yzone = mNode.maxY - mNode.minY;
	curPnt.clear();		curText.clear();
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
		
	    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    	{
	        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
	        if( geometry.valid() )
	        {
				osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
				for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
				{
					float u, v;
					u = ((pvx->x() - mNode.minX )/ Xzone) * WRAPTIMES;
					v = ((pvx->y() - mNode.minY )/ Yzone) * WRAPTIMES;
					TerPoint tp; tp.px = pvx->x(); tp.py = pvx->y(); tp.pz = pvx->z(); tp.nx = 0; tp.ny = 0; tp.nz = 128; curPnt.push_back(tp);
					TerTextcoord tt; tt.fU = u; tt.fV = v; curText.push_back(tt);
				}
			}//if valid
		}//n
	}//i
}


/* ��һ�� Group����ȹ�ߣ��������������ڵ����������֮ǰ��ҪԤ�ȵ��� CalculateMaxMinOfGroup() ������ mNode��*/
void CreateSkirtAndAddTextureCoordinateIntoSkirt(osg::ref_ptr<osg::Group> gp)
{
	/* �ȸ�һ���߻�����ȹ�ߡ�*/
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));

		/* Ȼ���ÿ���߻����ȹ�ߡ�*/
		for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
		{
			osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
			if( geometry.valid() )
			{
				for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
				{
					osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
					
					switch( pset->getType() )
					{
						case osg::PrimitiveSet::DrawArraysPrimitiveType:
						{
							switch(pset->getMode() )
							{
								case osg::PrimitiveSet::LINE_LOOP:
								{
									/* �� LINE_LOOP ���һ��ȹ�� */
									osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
									osg::ref_ptr<osg::Geometry> gm = ConvertFromLINELOOPToSkirtFace(gvx.get());
									if(gm != NULL) gnode->replaceDrawable( geometry, gm.get() );
								}
								break;
							}
						}// switch mode
						break;
					}//switch type
				}// k
			}//if valid
		}//n		
	}// i
	
	CreateVecArrayAndAddTextureCoordinateIntoFace(gp);
}

/* ��һ��Group ת���� DataStru�ṹ���Ա��������Ϊbtg�ļ� */
void CreateFaceDataStruFromGroup(osg::ref_ptr<osg::Group> gp)
{
	char sMater[64];

	/* ����Group���DataStru����������Ĳ��� */
	char cObjType;
	int indexType, indexOffset = 0;
	curTri.clear();
	cObjType = 10; indexType = 11;
	totalNormal=1;
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
		
		TriType tt;
		tt.mElement.clear();
		for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    {
			osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
			osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
			int vert_size = gvx->size();
			
			TriElement te; te.mVertex.clear();
			if( geometry.valid() )
			{
				for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
				{
					osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
			        
					switch( pset->getType() )
					{
						case osg::PrimitiveSet::DrawArraysPrimitiveType:
						{
							osg::DrawArrays *da = dynamic_cast<osg::DrawArrays *>(pset.get());
							int first = da->getFirst();
							unsigned int count = da->getCount();
							
							switch(pset->getMode() )
							{
								case osg::PrimitiveSet::TRIANGLES:
								{
									te.mVertex.clear();
									for(unsigned int ii=0; ii<count; ii++)
									{
										TriVertex tv;
										tv.P = indexOffset + ii; tv.N = 0; tv.T = indexOffset + ii; te.mVertex.push_back(tv);
									}
								}
								break;
							}
						}// switch mode
						break;
						case osg::PrimitiveSet::DrawElementsUIntPrimitiveType:
						{
							osg::DrawElementsUInt *deus = dynamic_cast<osg::DrawElementsUInt *>(pset.get());
							if( deus )
							{
								switch(pset->getMode() )
								{
									case osg::PrimitiveSet::TRIANGLES:
									{
										te.mVertex.clear();
										for( osg::DrawElementsUInt::iterator p = deus->begin(); p != deus->end(); p++ )
										{
											TriVertex tv;
											tv.P = indexOffset + *p; tv.N = 0; tv.T = indexOffset + *p; te.mVertex.push_back(tv);
										}
									}
									break;
								}
							}
						}
						break;
					}//switch type
				}// k
			}//if valid
			indexOffset += vert_size;
			te.iSum = te.mVertex.size();
			tt.mElement.push_back(te);
		}//n
		tt.iElems = tt.mElement.size();

		/* ��ȡ�������� */
		sprintf(sMater, "%s", (gnode->getName()).c_str());
		if(!strcmp(sMater, "IslandInOcean"))	sprintf(sMater, "Island");
		if((sMater[0] == 'L') || (sMater[0] == 'l'))		sprintf(tt.cMtype, "Lake");
		else if((sMater[0] == 'O') && (sMater[1] == 'c'))	sprintf(tt.cMtype, "Lake");
		else if((sMater[0] == 'R') || (sMater[0] == 'r'))	sprintf(tt.cMtype, "Road");
		else strcpy(tt.cMtype, sMater);

		tt.cObjType = cObjType;
		tt.cPropVal = indexType;
		curTri.push_back(tt);
	}//i
}

/* �����·��������ȹ��geode */
int CreateSkirtDataStruFromSkirtGroup(osg::ref_ptr<osg::Group> gp)
{
	/* ����gp���н�㣬�ҳ����ֵ����Сֵ */
	CalculateMaxMinOfGroup(gp);

	/* ��ȹ������������ꡣ*/
	CreateSkirtAndAddTextureCoordinateIntoSkirt(gp);
	
	/* ������DataStru�ṹ�� */
	CreateFaceDataStruFromGroup(gp);
	
	/* ��дbtg�ļ�ͷ */
	WriteBTGHeader();
	
	return 1;
}

/* �ռ���·���������������� */
void CreateVecArrayAndAddTextureCoordinateIntoRoad(osg::ref_ptr<osg::Group> gp)
{
	curPnt.clear();	curText.clear();
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
		
	    for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
    	{
	        osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
	        if( geometry.valid() )
	        {
				osg::ref_ptr<osg::Vec3Array> gvx   = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
				osg::ref_ptr<osg::Vec2Array> gtx   = dynamic_cast<osg::Vec2Array *>(geometry->getTexCoordArray(0));
				osg::Vec2Array::iterator t = gtx->begin();
				for( osg::Vec3Array::iterator pvx = gvx->begin(); pvx != gvx->end(); pvx++)
				{
					float u, v;
					if (t < gtx->end())	{	u = (*t)[0]; v = (*t)[1]; t++;	}

					/* ������Ϣ�����ݽṹ */
					TerPoint tp; tp.px = pvx->x(); tp.py = pvx->y(); tp.pz = pvx->z(); tp.nx = 0; tp.ny = 0; tp.nz = 128; curPnt.push_back(tp); 
					TerTextcoord tt; tt.fU = u; tt.fV = v; curText.push_back(tt);
				}
			}//if valid
		}//n
	}//i
}

/* ���������· */
int CreateRoadDataStruAndBtg(osg::ref_ptr<osg::Group> gp)
{
	/* ����gp���н�㣬�ҳ����ֵ����Сֵ */
	CalculateMaxMinOfGroup(gp);

	/* ��ʼ��ȫ������ */
	totalNormal = 0;

	/* ����·����������ꡣ*/
	CreateVecArrayAndAddTextureCoordinateIntoRoad(gp);
	
	/* ��··���ϵ�����Ƶ�Ҫ��������ɡ�ȡ���ڷ����������������������ƹ⡣2012-04-11 */
#if 0
	/* �ݲ����ɵ�·�ϵ�����Ƶ㡣��2013-01-22 */
	CreateCultPointsFromGroup(gp);
#endif
		
	/* ����Group���DataStru����������Ĳ��� */
	CreateFaceDataStruFromGroup(gp);

	/* дbtg�ļ�ͷ */
	WriteBTGHeader();
	
	return 1;
}

/* Add by Galen 2013-02-20 */
///////////////////////////// ���´�������� ////////////////////////
/* ���dclist */
void OutPutDclist(std::vector<osg::ref_ptr<ArealConstraint> > dclist)
{
	osg::ref_ptr<osg::Geode > gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
	{
		osgUtil::DelaunayConstraint * dc = dynamic_cast<osgUtil::DelaunayConstraint *> (p->get());
		if(dc->getNumPrimitiveSets() > 0)
		{
			osg::ref_ptr<osg::Geometry > gm = new osg::Geometry;
			gd->addDrawable(gm);
			gm->setVertexArray(dc->getVertexArray());
			gm->addPrimitiveSet(dc->getPrimitiveSet(0));

	    // Add an overall white color
	    osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
	    color->push_back(osg::Vec4(1,1,0,1));
	    gm->setColorArray(color.get());
	    gm->setColorBinding( osg::Geometry::BIND_OVERALL );
		}
	}
	osgDB::writeNodeFile(*gd, "F:\\dclist.ive");
}

/* ���dt */
void OutPutDt(osg::ref_ptr<osg::Vec3Array> coords, osg::ref_ptr<TerrainTriangulator> dt)
{
	osg::ref_ptr<osg::Geode > gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	// Create the new geometry from the retriangulated coordinates
	osg::ref_ptr<osg::Geometry > newGeometry = new osg::Geometry;
	newGeometry->setVertexArray( coords.get() );
	gd->addDrawable(newGeometry.get());
	
	// Add an overall white color
	osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
	color->push_back(osg::Vec4(1,1,0,1));
	newGeometry->setColorArray(color.get());
	newGeometry->setColorBinding( osg::Geometry::BIND_OVERALL );
	
	// Normals will be added by the smoothing visitor later
	newGeometry->addPrimitiveSet( dt->getTriangles() );
	osgDB::writeNodeFile(*gd, "F:\\dt.ive");
}

/* ���Points*/
void OutPutPoints(osg::ref_ptr<osg::Vec3Array> va)
{
	osg::ref_ptr<osg::Geode > gd = new osg::Geode;
	gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	
	for(osg::Vec3Array::iterator iter = va->begin(); iter != va->end(); iter++)
	{
		gd->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(iter->x(), iter->y(), iter->z()), 1.0)));
	}
	osgDB::writeNodeFile(*gd, "F:\\Point.ive");
}

///////////////////////////// ���ϴ�������� ////////////////////////


/* ����ͬ�����ཻ�ĵ������ */
#define UNITDISTANCE 300.0
void CreateNewStyleFaceOfLake(osg::ref_ptr<osg::Group> gp)
{
	osg::ref_ptr<osg::Group> newGp = new osg::Group;
	
	for(unsigned int m=0; m<gp->getNumChildren(); ++m)
	{
		osg::ref_ptr<osg::Geode> gFace = dynamic_cast<osg::Geode *>(gp->getChild(m));

#if 0
		osgDB::writeNodeFile(*gFace, "F:\\gFace.ive");
#endif

		for( unsigned int ii = 0; ii < gFace->getNumDrawables(); ii++ )
		{
			osg::ref_ptr<osg::Geometry> gm = dynamic_cast<osg::Geometry *>(gFace->getDrawable(ii));
			osg::ref_ptr<osg::Vec3Array> arpts = dynamic_cast<osg::Vec3Array *>(gm->getVertexArray());
			if(arpts->size() < 5) continue;

			/* �����߻��ļ��� */
			std::vector<osg::ref_ptr<osg::Vec3Array> > perList;
			TerrainProcessingUtils::computePerimeterAll( gm, perList);
			std::vector<osg::ref_ptr<ArealConstraint> > dclist;

			/* ���߻����dcѹ�뵽dclist */
			for( std::vector<osg::ref_ptr<osg::Vec3Array> >::iterator q= perList.begin(); q != perList.end(); q++ )
			{
				if (q->get()->size()>=3) 
				{
					osg::ref_ptr<ArealConstraint> dc = new ArealConstraint;
					dc->setVertexArray( q->get() );
					dc->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,q->get()->size()) );
					dclist.push_back( dc.get() );
				}
			}	//PerList
#if 0
		OutPutDclist(dclist);
#endif
	
			/* ���ȼ��㵱ǰ���� gnode �ľ�γ�ȼ�ֵ */
			MaxMinLonLat ml, *mmll = &ml;
			InitMMLLEx(mmll, "");
			double elevTotal = 0.0;
			int sum = 0;
	
			for( osg::Vec3Array::iterator r = arpts->begin(); r != arpts->end(); r++ )
			{
				osg::Vec3 pSpace = getSpaceVertex(*r);
				
				ProjectionMap pm;
				getProjectionMapByVertex(&pm, pSpace);
				SetMMLLEx(mmll, pm.lon, pm.lat);	elevTotal += pm.elev;	sum++;
			}
			double elevMid = elevTotal / sum;
	
			/* ����� 1 ����ռ�ݵľ�γ���Ƕ��� deltaLon, deltaLat */
			Geodesy p1, p2, geo;
			Carton car;
			p1._lon = mmll->minLong; p1._lat = mmll->maxLati; p1._elevation = 0;
			p2._lon = mmll->maxLong; p2._lat = mmll->maxLati; p2._elevation = 0;
			double distLon = distanceM(&p1, &p2);
			p1._lon = mmll->minLong; p1._lat = mmll->minLati; p1._elevation = 0;
			p2._lon = mmll->minLong; p2._lat = mmll->maxLati; p2._elevation = 0;
			double distLat = distanceM(&p1, &p2);
			double deltaLon = (mmll->maxLong - mmll->minLong) / distLon;
			double deltaLat = (mmll->maxLati - mmll->minLati) / distLat;
			
			/* ����� */
			for(double i = mmll->minLong; i < mmll->maxLong; i += UNITDISTANCE * deltaLon)
			{
				for(double j = mmll->minLati; j < mmll->maxLati; j += UNITDISTANCE * deltaLat)
				{
					/* ��γ��ת����xyz */
					geo._lon = i; geo._lat = j; geo._elevation = elevMid;	GeodToCart(&geo, &car);
					osg::Vec3 currVert(car.x, car.y, car.z);
		
					/* ����ÿ�����Ƿ�� g ���ཻ���󽻵� */
					osg::Vec3 isect;
					if( TerrainProcessingUtils::getIntersect( gFace, currVert, isect ))
					{
						car.x = isect.x();	car.y = isect.y(); car.z = isect.z();
						CartToGeod(&car, &geo);	geo._elevation = elevMid;
						GeodToCart(&geo, &car);
						osg::Vec3 cVert(car.x, car.y, car.z);
						arpts->push_back(cVert);
					}
				} // j
			} // i 

#if 0
			OutPutPoints(arpts);
#endif
			
			//�����������ǻ� ��������Լ��
			osg::ref_ptr<osg::Vec3Array> dcoords = TerrainProcessingUtils::uniqueCoords(arpts.get());		//( coords.get());
			osg::ref_ptr<TerrainTriangulator> dt = new TerrainTriangulator(dcoords.get());
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
				dt->addInputConstraint( p->get() );
			dt->triangulate();
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
				dt->removeInternalTriangles( p->get() );
			for( std::vector<osg::ref_ptr<ArealConstraint> >::iterator p= dclist.begin(); p != dclist.end(); p++ )
			{
				osg::ref_ptr<osg::Vec3Array> pts = p->get()->getPoints(dcoords.get());
				
				osg::ref_ptr<osg::Geode > gd = new osg::Geode;
				gd->setName("Lake");
				gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
				// Create the new geometry from the retriangulated coordinates
				osg::ref_ptr<osg::Geometry > newGeometry = new osg::Geometry;
				newGeometry->setVertexArray( pts );
				gd->addDrawable(newGeometry.get());
				
				// Normals will be added by the smoothing visitor later
				newGeometry->addPrimitiveSet( p->get()->getTriangles() );

#if 0
				// Add an overall white color
				osg::ref_ptr<osg::Vec4Array> color = new osg::Vec4Array;
				color->push_back(osg::Vec4(1,1,0,1));
				newGeometry->setColorArray(color.get());
				newGeometry->setColorBinding( osg::Geometry::BIND_OVERALL );

				osgDB::writeNodeFile(*gd, "F:\\dt.ive");
#endif
				newGp->addChild(gd);
			} // p
		}	// ii (gFace)
	} // m (gp)
	
	/* ��newGp ��� gp */
	gp->removeChildren(0, gp->getNumChildren());
	for(unsigned int m=0; m<newGp->getNumChildren(); ++m)
		gp->addChild(newGp->getChild(m));
	
	return;
}



/* ��Groupת����BTG�ļ������������͵ĺ����� */
int CreateNewStyleLakeDataStruAndBtg(osg::ref_ptr<osg::Group> gp)
{
	/* ���ԣ����������͵ĺ����� */
	CreateNewStyleFaceOfLake(gp);
	
	/* ����gp���н�㣬�ҳ����ֵ����Сֵ */
	CalculateMaxMinOfGroup(gp);

	/* ��ʼ��ȫ������ */
	totalNormal = 0;

	/* �ռ����㼯�Ϻ��������꼯�� */
	CreateVecArrayAndAddTextureCoordinateIntoFace(gp);

	/* ����Group���DataStru����������Ĳ��� */
	CreateFaceDataStruFromGroup(gp);

	/* дbtg�ļ�ͷ */
	WriteBTGHeader();
	
	return 1;
}
/* Add by Galen 2013-02-20 End */


/* ���������棬��Groupת����Btg�ļ��� */
int CreateLakeDataStruAndBtg(osg::ref_ptr<osg::Group> gp)
{
	/* ����gp���н�㣬�ҳ����ֵ����Сֵ */
	CalculateMaxMinOfGroup(gp);

	/* ��ʼ��ȫ������ */
	totalNormal = 0;

	/* �ռ����㼯�Ϻ��������꼯�� */
	CreateVecArrayAndAddTextureCoordinateIntoFace(gp);

	/* ����Group���DataStru����������Ĳ��� */
	CreateFaceDataStruFromGroup(gp);

	/* дbtg�ļ�ͷ */
	WriteBTGHeader();
	
	return 1;
}


/* ����һ��Group����DataStru���Ա��������Ϊbtg�ļ�������shp��ʽ�Ļ���Ϣ������ɭ�֡����Դ�ȵȵ����ɡ� */
int CreateDataStruFromGroup(osg::ref_ptr<osg::Group> gp)
{
	if(gp==NULL) return 0;

	//��ʼ���ڴ�
	curPnt.clear();		totalNormal = 0;
	curText.clear();	allLakeSkirt->removeChildren(0, allLakeSkirt->getNumChildren());
	curTri.clear();		allRoadSkirt->removeChildren(0, allRoadSkirt->getNumChildren());
	allDynamicGrass->removeChildren(0, allDynamicGrass->getNumChildren());
	allIslands->removeChildren(0, allIslands->getNumChildren());

	osg::ref_ptr<osg::Group> allForest = new osg::Group;
	osg::ref_ptr<osg::Group> allLake   = new osg::Group;
	osg::ref_ptr<osg::Group> allOcean  = new osg::Group;
	osg::ref_ptr<osg::Group> allLights = new osg::Group;
	osg::ref_ptr<osg::Group> allRoad   = new osg::Group;
	osg::ref_ptr<osg::Group> allRiver  = new osg::Group;
	osg::ref_ptr<osg::Group> allDetail = new osg::Group;
	osg::ref_ptr<osg::Group> allAirport= new osg::Group;
	
	/* ���㼫ֵ��mNode��Ϊ����ȡ��Χ�� */
	mNode.maxX = mNode.maxY = mNode.maxZ = -FLT_MAX;	mNode.minX = mNode.minY = mNode.minZ = FLT_MAX;
	
	/* �鿴group�����Ƿ���LINE_LOOP��������������ǻ� */
	bool isLineLoop = false;

	/* ���������е�Geode���ű��� */
	for(unsigned int i = 0; i < gp->getNumChildren(); i++)
	{
		osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));

		/* ������gnode�Ǻ���, ��������뵽allLake���С�*/
		/* ������gnode��ɭ��, ��������뵽allForest���С�*/
		/* ������gnode�ǵ�·, ��������뵽allRoad���С�*/
		/* ������gnode�Ǻ���, ��������뵽allRiver���С�*/
		if (gnode->getName().substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING) 
			allLake->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_OCEAN_STRING)) == DEFINE_OCEAN_STRING) 
			allOcean->addChild(gnode);
		if ((gnode->getName().substr(0, strlen(DEFINE_FOREST_STRING)) == DEFINE_FOREST_STRING) ||
			(gnode->getName().substr(0, strlen(DEFINE_FOREST_STRING2)) == DEFINE_FOREST_STRING2) )
			allForest->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)
			allRoad->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_RIVER_STRING)) == DEFINE_RIVER_STRING)
			allRiver->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_DETAIL_STRING)) == DEFINE_DETAIL_STRING)
			allDetail->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_LIGHTS_STRING)) == DEFINE_LIGHTS_STRING)
			allLights->addChild(gnode);
		if (gnode->getName().substr(0, strlen(DEFINE_AIRPORT_STRING)) == DEFINE_AIRPORT_STRING)
			allAirport->addChild(gnode);

		/* ������gnode�ǵ�·, ��������뵽allRoadSkirt���С�*/
		if (gnode->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING)
		{
			osg::ref_ptr<osg::Geode> gd = dynamic_cast<osg::Geode*>(gnode->clone(osg::CopyOp::DEEP_COPY_ALL));
			allRoadSkirt->addChild(gd);
		}

		for( unsigned int n = 0; n < gnode->getNumDrawables(); n++ )
		{
			osg::ref_ptr<osg::Geometry> geometry = dynamic_cast<osg::Geometry *>(gnode->getDrawable(n));
			if( geometry.valid() )
			{
				for(unsigned int k = 0; k< geometry->getNumPrimitiveSets(); k++)
				{
					osg::ref_ptr<osg::PrimitiveSet> pset = geometry->getPrimitiveSet(k);
					
					switch( pset->getType() )
					{
						case osg::PrimitiveSet::DrawArraysPrimitiveType:
						{
							switch(pset->getMode() )
							{
								case osg::PrimitiveSet::LINE_LOOP:
								{
									isLineLoop = true;
									/* �� LINE_LOOP ���ǻ���ʹ֮���һ���� */
									osg::ref_ptr<osg::Vec3Array> gvx = dynamic_cast<osg::Vec3Array *>(geometry->getVertexArray());
									osg::ref_ptr<osg::Geometry> gm = ConvertFromLINELOOPToTranglesFace(gvx.get());
									if(gm != NULL) gnode->replaceDrawable( geometry, gm.get() );
								}
								break;
							}
						}// switch mode
						break;
					}//switch type
				}// k
			}//if valid
		}//n

		/* ������gnode�Ǻ���, ��������뵽allLakeSkirt���С�*/
		if (gnode->getName() == DEFINE_LAKE_STRING)
		{
			osg::ref_ptr<osg::Geode> gd = dynamic_cast<osg::Geode*>(gnode->clone(osg::CopyOp::DEEP_COPY_ALL));
			allLakeSkirt->addChild(gd);
		}

		/* ������gnode�Ƕ�̬��, ��������뵽allDynamicGrass���С�*/
		if (gnode->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) == DEFINE_DYNAMICGRASS_STRING)
		{
			osg::ref_ptr<osg::Geode> gd = dynamic_cast<osg::Geode*>(gnode->clone(osg::CopyOp::DEEP_COPY_ALL));
			allDynamicGrass->addChild(gd);
		}

		/* ������gnode�Ƕ�̬��, ��������뵽allDynamicGrass���С�*/
		if (gnode->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING)
		{
			osg::ref_ptr<osg::Geode> gd = dynamic_cast<osg::Geode*>(gnode->clone(osg::CopyOp::DEEP_COPY_ALL));
			allIslands->addChild(gd);
		}
		
	}//i


#if 0			//�Ȳ���������ƽ������2012-07-20
	/* �����LINE_LOOP���������LINE_LOOP�Ǻ���,����һ��ƽ������ */
	if((isLineLoop)&&(allLake->getNumChildren() > 0))	
		AddCorrectionSpaceMesh(gp);
#endif

	/* �����ɭ�֣�����Ҫ��������*/
	if(allForest->getNumChildren() > 0)
	{
		CreateTreePointsFromGroup(allForest);
		
		/* ��gp���������ɭ�ֵ�geode�Ƴ� */
		osg::ref_ptr<osg::Group> gpNoForest = new osg::Group;
		for(unsigned int i = 0; i < gp->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
			if (!((gnode->getName().substr(0, strlen(DEFINE_FOREST_STRING)) == DEFINE_FOREST_STRING) || (gnode->getName().substr(0, strlen(DEFINE_FOREST_STRING2)) == DEFINE_FOREST_STRING2)))
				gpNoForest->addChild(gnode);
		}
		if(gpNoForest->getNumChildren() == 0) return 0;

		gp->removeChildren(0, gp->getNumChildren());
		for(unsigned int i = 0; i < gpNoForest->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gpNoForest->getChild(i));
			gp->addChild(gnode);
		}		
	}

	/* �������ƣ�Ҳ��Ҫ��������*/
	if(allLights->getNumChildren() > 0)
	{
		CreateCultPointsFromGroup(allLights);
		
		/* ��gp�����������Ƶ�geode�Ƴ� */
		osg::ref_ptr<osg::Group> gpNoLights = new osg::Group;
		for(unsigned int i = 0; i < gp->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
			if (!(gnode->getName().substr(0, strlen(DEFINE_LIGHTS_STRING)) == DEFINE_LIGHTS_STRING))
				gpNoLights->addChild(gnode);
		}
		if(gpNoLights->getNumChildren() == 0) return 0;

		gp->removeChildren(0, gp->getNumChildren());
		for(unsigned int i = 0; i < gpNoLights->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gpNoLights->getChild(i));
			gp->addChild(gnode);
		}		
	}

	/* ����ǵ�·����ô�����������㡣����Ǳ���Ļ���Ϣ��������������һ���� */
	if(allRoad->getNumChildren() > 0)
	{
		/* ��gp��������е�·��geode�Ƴ� */
		osg::ref_ptr<osg::Group> gpNoRoads = new osg::Group;
		for(unsigned int i = 0; i < gp->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
			if (!(gnode->getName().substr(0, strlen(DEFINE_ROAD_STRING)) == DEFINE_ROAD_STRING))
				gpNoRoads->addChild(gnode);
		}
		if(gpNoRoads->getNumChildren() == 0) return 0;

		gp->removeChildren(0, gp->getNumChildren());
		for(unsigned int i = 0; i < gpNoRoads->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gpNoRoads->getChild(i));
			gp->addChild(gnode);
		}		
	}

	/* ����Ǻ�����ôBTG�ļ����㡣����Ǳ���Ļ���Ϣ��������������һ���� */
	/* ˳��Ѷ�̬��Ҳһ��ȥ������*/
	if((allLake->getNumChildren() > 0) ||
		 (allDynamicGrass->getNumChildren() > 0) ||
		 (allIslands->getNumChildren() > 0) )
	{
		/* ��gp��������к�����geode�Ƴ� */
		osg::ref_ptr<osg::Group> gpNoLakes = new osg::Group;
		for(unsigned int i = 0; i < gp->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gp->getChild(i));
			if((!(gnode->getName().substr(0, strlen(DEFINE_LAKE_STRING)) == DEFINE_LAKE_STRING)) &&
				 (!(gnode->getName().substr(0, strlen(DEFINE_DYNAMICGRASS_STRING)) == DEFINE_DYNAMICGRASS_STRING)) &&
				 (!(gnode->getName().substr(0, strlen(DEFINE_ISLAND_STRING)) == DEFINE_ISLAND_STRING)) )
				gpNoLakes->addChild(gnode);
		}
		if(gpNoLakes->getNumChildren() == 0) return 0;

		gp->removeChildren(0, gp->getNumChildren());
		for(unsigned int i = 0; i < gpNoLakes->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gpNoLakes->getChild(i));
			gp->addChild(gnode);
		}	
	}
	
	/* �������н�㣬�ҳ����ֵ����Сֵ */
	CalculateMaxMinOfGroup(gp);

	/* �ռ����㼯�Ϻ��������꼯�� */
	CreateVecArrayAndAddTextureCoordinateIntoFace(gp);

	/* ����Group���DataStru����������Ĳ��� */
	CreateFaceDataStruFromGroup(gp);

	/* дbtg�ļ�ͷ */
	WriteBTGHeader();
	
	return 1;
}


/* �ж�һ����������˳ʱ�뻹����ʱ�� 2012-10-18 */
/* ��ʱ�뷵����ֵ��˳ʱ�뷵�ظ�ֵ�����㹲�߷���0 */
double CounterClockWise(osg::Vec3 a, osg::Vec3 b, osg::Vec3 c)
{
	return ((b.x() - a.x())*(c.y() - b.y()) - (c.x() - b.x())*(b.y() - a.y()));   
}

/*******************************************************
 *
 *  ����һ���߶����һ�����Σ�����shp��ʽ�Ļ���Ϣ�ĵ�·�����ɡ�
 *
 *******************************************************/
/* ��ȡһ���߶ε�ƽ�����ϵ����������㣬����ԭ�е����������һ�����Σ����������Σ�������ɲ������� */
/* 2012-10-18 ԭ�����㷨��ֱ�Ӹ��ݿռ��ϵ����������߶�����չ�ɾ��Σ��������������˾��ββ�룬�����ܻᷢ������*/
/* ����֮����㷨�����ռ����ֱ������ϵת��Ϊ��γ������ϵ�������߶�ͳһΪ���θ߶�0.0�������ƽ���ڣ�ʹ�þ�γ������ϵ��ֵͳһ���㡣*/
osg::ref_ptr<osg::Geometry> CreateRectangleWithEdge(SimEdge *se, float Width, int Direction)
{
	osg::Vec3 V[2];
	double xWidth, yWidth;

	/* �Ƚ�ֱ������ϵת��Ϊ��γ������ϵ */
	Geodesy gA, gB, gV0, gV1;
	Carton cA, cB, cV0, cV1;
	cA.x = se->A.x();	cA.y = se->A.y(); cA.z = se->A.z();	CartToGeod(&cA, &gA);
	cB.x = se->B.x();	cB.y = se->B.y(); cB.z = se->B.z();	CartToGeod(&cB, &gB);

	/* ����� 1 ����ռ�ݵľ�γ���Ƕ��� deltaLon, deltaLat */
	Geodesy p1, p2;
	p1._lon = (int)(gA._lon);		p1._lat = gA._lat;	p1._elevation = 0;
	p2._lon = (int)(gA._lon + 1);	p2._lat = gA._lat;	p2._elevation = 0;
	double distLon = distanceM(&p1, &p2);
	p1._lon = gA._lon; p1._lat = (int)(gA._lat);		p1._elevation = 0;
	p2._lon = gA._lon; p2._lat = (int)(gA._lat + 1);	p2._elevation = 0;
	double distLat = distanceM(&p1, &p2);
	double deltaLon = 1.0 / distLon;
	double deltaLat = 1.0 / distLat;

	/* �������֪�߶ε�б�� */
	if(gA._lon != gB._lon)
	{
		double yDelta = gA._lat - gB._lat;
		double xDelta = gA._lon - gB._lon;
		double k = yDelta / xDelta;
		
		/* ֱ����Ƕ� */
		double alpha  = atan(k);
		double alphaC = alpha * 180.0 / osg::PI;
		
		alpha += 0.5 * osg::PI;
		double alphaC2 = alpha * 180.0 / osg::PI;

		/* ����б�ʺͿ�ȼ���x,y�����ƫ�ƣ�����Ƕȣ�Ȼ����sin,cos */
		yWidth = Width * deltaLat * sin(alpha);
		xWidth = Width * deltaLon * cos(alpha);
		if(xDelta < 0.0)
		{
			yWidth = - yWidth; xWidth = - xWidth;
		}
	}else{
		xWidth = Width * deltaLon; yWidth = 0;
	}

	if(Direction)	{	xWidth = - xWidth; yWidth = - yWidth;	}

	/* ���������˵������������������ */
	gV0._lon =  gA._lon + xWidth;	gV0._lat = gA._lat  + yWidth;	gV0._elevation = 0.0;	GeodToCart(&gV0, &cV0);
	gV1._lon =  gB._lon + xWidth;	gV1._lat = gB._lat  + yWidth;	gV1._elevation = 0.0;	GeodToCart(&gV1, &cV1);
	V[0].x() = cV0.x;	V[0].y() = cV0.y;	V[0].z() = cV0.z;
	V[1].x() = cV1.x;	V[1].y() = cV1.y;	V[1].z() = cV1.z;
	
	/* �ж�һ����������˳ʱ�뻹����ʱ�� 2012-10-18 */
	double clockwire = CounterClockWise(se->A, se->B, V[0]);
	if((clockwire > 0.0)&&(Direction == 0) || (clockwire < 0.0)&&(Direction == 1))
	{
		gV0._lon =  gA._lon - xWidth;	gV0._lat = gA._lat  - yWidth;	gV0._elevation = 0.0;	GeodToCart(&gV0, &cV0);
		gV1._lon =  gB._lon - xWidth;	gV1._lat = gB._lat  - yWidth;	gV1._elevation = 0.0;	GeodToCart(&gV1, &cV1);
		V[0].x() = cV0.x;	V[0].y() = cV0.y;	V[0].z() = cV0.z;
		V[1].x() = cV1.x;	V[1].y() = cV1.y;	V[1].z() = cV1.z;
	}

	/* ����һ��Geometry*/	
	osg::ref_ptr<osg::Geometry> gm = new osg::Geometry;
	osg::ref_ptr<osg::Vec3Array> vertex = new osg::Vec3Array;
	vertex->push_back(se->A);
	vertex->push_back(se->B);
	vertex->push_back(V[0]);
	vertex->push_back(V[1]);
	gm->setVertexArray(vertex);
	osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	dui->push_back(0); dui->push_back(1); dui->push_back(2);
	dui->push_back(3); dui->push_back(2); dui->push_back(1);
	gm->addPrimitiveSet(dui.get());
#if 0		//Use color
	osg::Vec4Array* colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f,1.0f,0.0f,1.0f));
    gm->setColorArray(colors);
    gm->setColorBinding(osg::Geometry::BIND_OVERALL);
#endif
	return gm.release();	
}

/*******************************************************
 *
 *  ��ȡһ�鶥��ľ�γ�����꼫ֵ
 *
 *******************************************************/

void InitMMLLEx(MaxMinLonLat *m, string fn)
{
	m->mShpName = fn;
	m->minLong = m->minLati = 180.0;
	m->maxLong = m->maxLati = -180.0;
}

void SetMMLLEx(MaxMinLonLat *m, double lon, double lat)
{
	if(lon > m->maxLong) m->maxLong = lon;
	if(lat > m->maxLati) m->maxLati = lat;
	if(lon < m->minLong) m->minLong = lon;
	if(lat < m->minLati) m->minLati = lat;
}


/*******************************************************
 *
 *  ��ȡһ���Ƶ��ļ����������
 *
 *******************************************************/
std::vector<LightNode> allLights;
std::vector<LightNode> directionalLights;
std::vector<LightNode> undirectionalLights;
std::vector<string>	   allNames;
std::vector<ColorName> allColors;
bool foundLight = false;

/*******************************************************
 *
 *  ���¹����Ƕ�ȡ�û�flt�����ĵƵ㡣
 *
 *******************************************************/

/* ����flt�ļ��������еĵƵ㣬����Щ�Ƶ㱣���� allLights���档*/
class FindLightVisitor : public osg::NodeVisitor
{
public:
	FindLightVisitor() : osg::NodeVisitor(TRAVERSE_ALL_CHILDREN) { }; 

	void apply(osg::Group &gp)
	{
		if (gp.getName().substr(0, strlen("Lights")) == "Lights")
		{
			foundLight = true;
			doApply(gp);
		}
		traverse(gp);
	}

	void readLightPointNode(osg::ref_ptr<osgSim::LightPointNode> lpn, LightNode *ln, osg::Matrix m)
	{
		if(lpn)
		{
			LightNode lnd;
			if(ln == NULL)
			{
				lnd.name = lpn->getName();	lnd.vlist.clear();	allNames.push_back(lnd.name);
				allLights.push_back(lnd);
				ln = &lnd;
			}

			for(unsigned int j = 0; j<lpn->getNumLightPoints(); j++)
			{
				LightVertex lv;
				osgSim::LightPoint lp =  lpn->getLightPoint(j);
				lv.position = lp._position * m;
				lv.color    = lp._color;
				
				/* ��ȡ�Ƶ�ķ���2012-03-26 */
				osgSim::Sector * sec = (osgSim::Sector *)(lp._sector);
				osg::ref_ptr<osgSim::DirectionalSector> ds   = dynamic_cast<osgSim::DirectionalSector *>(sec);	//reinterpret_cast
				if(ds) lv.direction = ds->getDirection(); else lv.direction.set(0, 0, 1);
				ln->vlist.push_back(lv);
				
				/* ����һ���µ���ɫ */
				int iFlag = 0;
				for(unsigned int k = 0; k<allColors.size(); k++)	{	if((lv.color == allColors[k].color)&&(ln->name == allColors[k].name)) { iFlag = 1; break; }}
				if(!iFlag) {	ColorName cn; cn.color = lv.color;	cn.name = ln->name;	allColors.push_back(cn);	}
			}
			unique = *ln;
		}
	}
	
	bool readMatrixTransform(osg::MatrixTransform* mt, LightNode *ln)
	{
		if(mt)
		{
			LightNode lnd;
			if(ln == NULL)
			{
				lnd.name = mt->getName();	lnd.vlist.clear();	allNames.push_back(lnd.name);
				allLights.push_back(lnd);
				ln = &lnd;
			}
			
			osg::Matrix m = mt->getMatrix();
			for(unsigned int k = 0; k<mt->getNumChildren(); k++)
			{
				osg::ref_ptr<osgSim::LightPointNode> lpn = dynamic_cast<osgSim::LightPointNode *>(mt->getChild(k));
				readLightPointNode(lpn, ln, m);
			}
			return true;
		}else
			return false;
	}

	bool readGroup(osg::Group* gp)
	{
		if(gp)
		{
			LightNode ln;
			ln.name = gp->getName();	ln.vlist.clear();	allNames.push_back(ln.name);
			for(unsigned int k = 0; k<gp->getNumChildren(); k++)
			{
				osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(gp->getChild(k));
				if(readMatrixTransform(mt, &ln) == false)
				{
					osg::ref_ptr<osgSim::LightPointNode> lpn = dynamic_cast<osgSim::LightPointNode *>(gp->getChild(k));
					osg::Matrix *mx = new osg::Matrix;
					readLightPointNode(lpn, &ln, *mx);
				}
			}
			allLights.push_back(ln);
			return true;
		}else
			return false;
	}

	void doApply(osg::Group &gp)
	{
		allLights.clear();	allNames.clear();	allColors.clear();
	
		for(unsigned int i = 0; i<gp.getNumChildren(); i++)
		{
			osg::Group* grp = dynamic_cast<osg::Group *>(gp.getChild(i));
			if(readGroup(grp) == false)
			{
				osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(gp.getChild(i));
				if(readMatrixTransform(mt, NULL) == false)
				{
					osg::ref_ptr<osgSim::LightPointNode> lpn = dynamic_cast<osgSim::LightPointNode *>(gp.getChild(i));
					osg::Matrix *mx = new osg::Matrix;
					readLightPointNode(lpn, NULL, *mx);
				}
			}
		}
	}

private:
	LightNode unique;	
};


/* ��matname���»�����Ϊ�ָ������ֳ�����С�Σ�������materialParts���� */
int splitMaterialName(string matname)
{
	char matpart[128], mat[256];
	int j = 0;
	
	/* ��matname���»�����Ϊ�ָ������ֳ�����С�Σ�������materialParts���� */
	sprintf(mat, "%s", matname.c_str());
	materialParts.clear();
	for(unsigned int i = 0; i<strlen(mat); i++)
	{
		if((mat[i] != '_')&&(mat[i] != '-')) { matpart[j++] = mat[i];	}
		if((mat[i] == '_')||(mat[i] == '-')) { matpart[j] = 0;	materialParts.push_back(matpart);	j=0; }
	}
	matpart[j] = 0;	materialParts.push_back(matpart);	j=0;
	
	return materialParts.size();
}

/* ��matname��б�ܻ��߷�б����Ϊ�ָ������ֳ�����С�Σ�������materialParts���� */
int splitPathName(std::string matname)
{
   char matpart[128], mat[256];
   int j = 0;
   
   /* ��matname��б�ܻ��߷�б����Ϊ�ָ������ֳ�����С�Σ�������materialParts���� */
   sprintf(mat, "%s", matname.c_str());
   materialParts.clear();
   for(unsigned int i = 0; i<strlen(mat); i++)
   {
   	if((mat[i] != '/')&&(mat[i] != '\\')) { matpart[j++] = mat[i];	}
   	if((mat[i] == '/')||(mat[i] == '\\')) { matpart[j] = 0;	materialParts.push_back(matpart);	j=0; }
   }
   matpart[j] = 0;	materialParts.push_back(matpart);	j=0;
   
   return materialParts.size();
}

/* �ϲ�matname,��materialParts�����ΰ�˳������� */
std::string mergeMaterialName()
{
	std::string matn; matn = materialParts[0];
	for(unsigned int i=1; i<materialParts.size(); i++)
	{
		matn = matn + "_" + materialParts[i];
	}
	return matn;
}

/* Ѱ�Ҷ���ƣ�2012-04-14 ������ */
/* ��AllLights���� �ռ����еĶ���ƣ����Ұ��վ���Զ������ */
/* �ٶ����ж���ƶ��ǳ��ŵģ�����������֮��ľ������10�ף������һ�ŵ���һ�顣*/

std::vector<LightNode> allDLGroups;
#define VHEIGHT 10.0

int GetDirectionalLights()
{
	unsigned int i, j, k, m;
	osg::Vec3d noDirection; 	noDirection.set(0.0, 0.0, 1.0);

	/* �����ռ����еĶ���ƣ�������directionalLights����*/
	/* �Ȱ�allLights����ΪdirectionalLights��undirectionalLights */
	directionalLights.clear();
	undirectionalLights.clear();
	for(i=0; i<allLights.size(); i++)
	{
		LightNode lnd;
		lnd.name = allLights[i].name;	lnd.vlist.clear();

		for(j=0; j<allLights[i].vlist.size(); j++)
		{
			LightVertex lv;
			lv.position  = allLights[i].vlist[j].position;
			lv.direction = allLights[i].vlist[j].direction;
			lv.color     = allLights[i].vlist[j].color;
			lv.index	 = -1;
			lnd.vlist.push_back(lv);
		}

		if(allLights[i].vlist[0].direction != noDirection)
		{
			directionalLights.push_back(lnd);
		}else{
			undirectionalLights.push_back(lnd);
		}
	}
	
	/* �ٰ�undirectionalLightsд��allLights */
	allLights.clear();	allColors.clear();
	for(i=0; i<undirectionalLights.size(); i++)
	{
		LightNode lnd;
		lnd.name = undirectionalLights[i].name;	lnd.vlist.clear();

		for(j=0; j<undirectionalLights[i].vlist.size(); j++)
		{
			LightVertex lv;
			lv.position  = undirectionalLights[i].vlist[j].position;
			lv.direction = undirectionalLights[i].vlist[j].direction;
			lv.color     = undirectionalLights[i].vlist[j].color;
			lv.index	 = -1;
			lnd.vlist.push_back(lv);
		}
		allLights.push_back(lnd);
		ColorName cn; cn.color = lnd.vlist[0].color;	cn.name = lnd.name;	allColors.push_back(cn);
	}
	
	/* ����Щ����ư������ƣ�ͬ���Ƶİ��վ�����顣һ������һ�Ŷ���ƶ��ǳ�һ������ġ�ͬ��֮�ڵĶ������̾����ǱȽϽ��� */
	/* ����Լ�������һ����ͬһ���֮�����̾���С��10�ף���ô����ƾ���������ơ�*/
	for(i=0; i<directionalLights.size(); i++)
	{
		j = 0;
		while(1)
		{
			/* ���ȣ����б����ҳ���һ��û����� LightVertex */
			while(j < directionalLights[i].vlist.size()) { if(directionalLights[i].vlist[j].index != -1) j++; else break; }
			if(j >= directionalLights[i].vlist.size()) break;
			
			/* �����LightVertex ���뵽 LightNode ���棬���ı�index��״̬��ֵ */
			LightNode dlg; dlg.vlist.clear(); dlg.name = directionalLights[i].name;
			LightVertex lv;
			lv.position  = directionalLights[i].vlist[j].position;
			lv.direction = directionalLights[i].vlist[j].direction;
			lv.color     = directionalLights[i].vlist[j].color;
			lv.index	 = 0;
			dlg.vlist.push_back(lv);
			directionalLights[i].vlist[j].index = 1;
			
			/* ��Ѱʣ��δ����ĵƵ㣬���ͬdlg�����еƵ����̾���С��10.0���Ͱ�����Ƶ����dlg���ڡ�*/
			for(k = j+1; k < directionalLights[i].vlist.size(); k++)
			{
				if(directionalLights[i].vlist[k].index == -1)
				{
					/* ������δ�����k����dlg�ڸ�����֮�����С���룬������min_length���档*/
					double min_length = 99999999.9;
					for(m=0; m<dlg.vlist.size(); m++)
					{
						osg::Vec3 vdist = dlg.vlist[m].position - directionalLights[i].vlist[k].position;
						double len = length(vdist);
						if(len < min_length) min_length = len;
					}
					
					/* �ж� */
					if(min_length < 10.0)
					{
						LightVertex lv;
						lv.position  = directionalLights[i].vlist[k].position;
						lv.direction = directionalLights[i].vlist[k].direction;
						lv.color     = directionalLights[i].vlist[k].color;
						lv.index	 = 0;
						dlg.vlist.push_back(lv);
						directionalLights[i].vlist[k].index = 1;
					}					
				}
			}	//for(k)
			
			/* ������ֺõ�����뵽allDLGroups���� */
			allDLGroups.push_back(dlg);
		}	//while (j)
	}
	
	/* ����allDLGroups�����ÿ���飬��������ǵķ���*/
	/* ��ÿ�����ڣ�ȡ�������㣬����������������߶εĴ�ֱƽ���ߣ������߶�����Ĵ�ֱƽ������ȡ����V[0], V[1] */
	/* �������֪�߶ε��е� */
	osg::Vec3d V[2];
	std::vector<string>	 singleNames;	singleNames.clear();
	for(i=0; i<allDLGroups.size(); i++)
	{
#if 0	/* ԭ�������������Ǹ���һ�Ŷ�����������������ƵĴ�ֱƽ������������ġ����������б䣺����ÿһ������ƶ����Լ��̶��ķ���
           ����ԭ���������Ͳ��ʵ��ˡ�
         */
         
         
		/* ������ڵƹ���ֻ��һ����û������������ڶ���Ƶķ���*/
		if(allDLGroups[i].vlist.size() < 2) continue;
			
		SimEdge lse; lse.A = allDLGroups[i].vlist[0].position; lse.B = allDLGroups[i].vlist[1].position; 
		SimEdge *se; se = &lse;
		osg::Vec3d mid = (se->A + se->B) / 2;
		
		/* �������֪�߶ε�б�ʣ�ֻ��ƽ���ڵ�б�� */
		if(se->A[0] != se->B[0])
		{
			float k = (se->A[1] - se->B[1]) / (se->A[0] - se->B[0]);
			
			/* ������ֱ�ֱཻ�ߵ�б����˻�Ϊ-1��k1*k2=-1. */
			float kv;
			if( k != 0.0) 	kv = -(1.0 / k);
			
			/* �����ֱƽ������������V[0], V[1] */
			V[0].x() = mid.x() + VHEIGHT; V[1].x() = mid.x() - VHEIGHT;
			if( k != 0.0)
			{
				V[0].y() = kv * (V[0].x() - mid.x()) + mid.y(); 
				V[1].y() = kv * (V[1].x() - mid.x()) + mid.y(); 
			}else{
				V[0].x() = V[1].x() = mid.x();
				V[0].y() = mid.y() + VHEIGHT; V[1].y() = mid.y() - VHEIGHT;
			}
			V[0].z() = mid.z(); V[1].z() = mid.z();
		}else{
			V[0].x() = mid.x() + VHEIGHT; V[1].x() = mid.x() - VHEIGHT;
			V[0].y() = V[1].y() = mid.y();
			V[0].z() = mid.z(); V[1].z() = mid.z();
		}
	
		/* ���ݴ�ֱƽ�����ϵ������������㷽�� */
		osg::Vec3 cart1, cart2;
		if(allDLGroups[i].vlist[0].direction.y() > 0.0)
		{
			cart1 = V[0];	cart2 = V[1];
		}else{
			cart1 = V[1];	cart2 = V[0];
		}
	    
	    double mLength;
	    float  angle = 3.0;
		
	    osg::Vec3 up = normalize(cart1);
	    // angle up specified amount
	    osg::Vec3 rwy_vec = cart2 - cart1;
	    mLength = length(rwy_vec);
	    double up_length = mLength * tan(angle * SGD_DEGREES_TO_RADIANS);
	    osg::Vec3 light_vec = normalize(rwy_vec + (up * up_length));
	    
	    /* ���������еĵƵ㷽������Ϊ light_vec */
		for(j=0; j<allDLGroups[i].vlist.size(); j++)
		{
			allDLGroups[i].vlist[j].direction = light_vec * rotateMatrix;
		}

#else	/* ���ڵ�������ֱ�Ӹ���ÿ������Ƶ�directionȡ����ֵ�������������Ƶķ���ǰ���Ǹ�������Ƶķ����Ƕ����ġ�*/
		for(j=0; j<allDLGroups[i].vlist.size(); j++)
		{
		    osg::Vec3 light_vec = normalize(allDLGroups[i].vlist[j].direction);
			allDLGroups[i].vlist[j].direction = normalize(light_vec * rotateMatrix);
		}
#endif

		
		/* �޸�allDLGGroups[i]�������������޸�Ϊ�� ??_??_DIRECTIONAL_??_??��  */
		splitMaterialName(allDLGroups[i].name);
		materialParts[2] = "DIRECTIONAL";
		allDLGroups[i].name = mergeMaterialName();
		
		/* ��allDLGGroups[i].name �������� singleNames, ͬ�������ֽ���һ�Ρ� */
		int iFlag = 0;
		for(j=0; j<singleNames.size(); j++)
		{
			if(singleNames[j] == allDLGroups[i].name)	{ iFlag = 1; break; }
		}
		if(!iFlag) singleNames.push_back(allDLGroups[i].name);
	}

    /* ��allDLGGroups�����ֵд�ص�allLights����ͬ������Ҫ������Ϊһ�� */
    for(i=0; i<singleNames.size(); i++)
    {
		LightNode lnd;
		lnd.name = singleNames[i];	lnd.vlist.clear();

		for(j=0; j<allDLGroups.size(); j++)
		{
			if(allDLGroups[j].name == singleNames[i])
			{
				for(k=0; k<allDLGroups[j].vlist.size(); k++)
				{
					LightVertex lv;
					lv.position  = allDLGroups[j].vlist[k].position;
					lv.direction = allDLGroups[j].vlist[k].direction;
					lv.color     = allDLGroups[j].vlist[k].color;
					lv.index	 = -1;
					lnd.vlist.push_back(lv);
				}
			}
		}

		allLights.push_back(lnd);
		ColorName cn; cn.color = lnd.vlist[0].color;	cn.name = lnd.name;	allColors.push_back(cn);
    }

	return 1;
}




/* ��ȡ����Flt�ļ��������������еĵƵ���Ϣ����д��btg���ݽṹ���档*/
int ReadFltAndCreateLightsDataStru(char *fname)
{
	unsigned int i, j, n;
	char apt[256];

	/* ��ȡǰ�����ɵĵƹ���ʱ�ļ���AptLight.ive�� */
	sprintf(apt, "%s\\AptLight.ive", tempDirectory);

	/* �������еĵƵ� */
	osg::ref_ptr<osg::Node> ap = osgDB::readNodeFile(apt);
	if(!ap.valid()) return 0;
	
	/* ɾ��ǰ�����ɵĵƹ���ʱ�ļ���AptLight.ive�� */
	DeleteFile(apt);

	FindLightVisitor flv;
	ap->accept(flv);

	/* ����Ҳ����ƹ⣬����0 */
	if(foundLight == false) return 0;
		
	/* ���㶨��� */
	GetDirectionalLights();
	
	/* ���ݲ�ͬ����ɫ���ҵ���ͬ�Ĳ������ƣ�дDataStru */
	/* �Ƚ�������ת�� */
	mNode.maxX = mNode.maxY = mNode.maxZ = -FLT_MAX;	mNode.minX = mNode.minY = mNode.minZ = FLT_MAX;
	int Idx = 0;
	curPnt.clear(); totalNormal = 0;
	for(i=0; i<allLights.size(); i++)
	{
		for(j=0; j<allLights[i].vlist.size(); j++)
		{
			osg::Vec3d v; v = allLights[i].vlist[j].position * coordMatrix;
			allLights[i].vlist[j].position.set(v);	allLights[i].vlist[j].index = Idx++;

			/* ����ÿ����ķ��������������ķ����� 2012-03-26 */
			TerPoint tp; tp.px = v.x(); tp.py = v.y(); tp.pz = v.z(); 
			osg::Vec3d dir; dir = allLights[i].vlist[j].direction; dir.normalize();
			tp.nx = (short)(dir.x() * 128); tp.ny = (short)(dir.y() * 128); tp.nz = (short)(dir.z() * 128);
			totalNormal++;
			curPnt.push_back(tp);

			/*��ֵ������mNode ����*/
			if(mNode.maxX < v.x()) mNode.maxX = v.x();
			if(mNode.maxY < v.y()) mNode.maxY = v.y();
			if(mNode.maxZ < v.z()) mNode.maxZ = v.z();
			if(mNode.minX > v.x()) mNode.minX = v.x();
			if(mNode.minY > v.y()) mNode.minY = v.y();
			if(mNode.minZ > v.z()) mNode.minZ = v.z();
		}
	}

	/* ͳ�Ʋ������� */
	for(i=0; i<allColors.size(); i++)
	{
		if(allColors[i].name.length() < 15 )
		{
			int idx = FindIdxByColor(allColors[i].color);
			if(idx != -1)	
			{
				for(j = 0; j < allLights.size(); j++)
				{
					if(allLights[j].name == allColors[i].name)
						allLights[j].name = colors[idx].name;
				}
				allColors[i].name = colors[idx].name;
			}
		}
	}
	
	/* ����������Ψһ�� */
	std::vector<ColorName> uniqueColor; uniqueColor.clear();
	for(i=0; i<allColors.size(); i++)
	{
		int iFlag = 0;
		for(j=0; j<uniqueColor.size(); j++)
		{
			if((allColors[i].name == uniqueColor[j].name) && (allColors[i].color == uniqueColor[j].color))
			{ iFlag = 1; break;	}
		}
		if(!iFlag) {	uniqueColor.push_back(allColors[i]);	}
	}	
	
	/* д����� */
	char cObjType;
	int indexType, indexOffset = 0;
	curTri.clear();
	cObjType = 9; indexType = 3;

	for(n = 0; n < uniqueColor.size(); n++)
	{
		TriType tt;
		tt.mElement.clear();

		/* ����te */
		TriElement te; te.mVertex.clear();
		for(i=0; i<allLights.size(); i++)
		{
			for(j=0; j<allLights[i].vlist.size(); j++)
			{
    		if((allLights[i].vlist[j].color == uniqueColor[n].color) && (allLights[i].name == uniqueColor[n].name))
    		{
					TriVertex tv;
					tv.P = tv.N = tv.T = allLights[i].vlist[j].index;  te.mVertex.push_back(tv);
    		}
			}
		}
		te.iSum = te.mVertex.size();
		tt.mElement.push_back(te);
		tt.iElems = tt.mElement.size();

		/* д���� */
		sprintf(tt.cMtype, "%s", uniqueColor[n].name.c_str());
		tt.cObjType = cObjType;
		tt.cPropVal = indexType;
		curTri.push_back(tt);
	}
	
	/* ����BTG�ļ�ͷ */
	WriteBTGHeader();

	return 1;
}


