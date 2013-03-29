#include "stdafx.h"
#include "d3dtypes.h"
#include "global.h"
#include "ShpFile.h"

#define USE_INSERT_VERTEX 1
#define INSERTDISTANCE 0.002

//外部函数调用及全局变量
extern float  mRoadWidthInShp;
float  mLightPotDistance;
extern double GetElevation(double lon, double lat, double pre);
extern void CalculateMaxMinOfGroup(osg::ref_ptr<osg::Group> gp);
extern int CreateDataStruFromGroup(osg::ref_ptr<osg::Group> gp);
extern osg::ref_ptr<osg::Geometry> CreateRectangleWithEdge(SimEdge *se, float Width, int Direction);
extern void intersect(  osg::Vec3 A, osg::Vec3 B, osg::Vec3 C, osg::Vec3 D, osg::Vec3 &I );
extern void interpolation(osg::Vec3 A, osg::Vec3 B, std::vector<osg::Vec3> &vec_list, float segs);	
extern int ifHasIslandInCurrentTerrain;
extern std::vector<std::string> allDetailedTxt;
extern char tempDirectory[_MAX_PATH];
extern char debugFile[_MAX_PATH];
extern bool CurrentLampType;
extern string getTreeName(string curr);
extern double LenOf2V(const osg::Vec3 v0, const osg::Vec3 v1);

//本地置顶变量及函数
MaxMinLonLat currLayer;
void SetMMLL(double lon, double lat);
typedef std::vector<SimEdge> MyPolyLine;
std::vector<MyPolyLine> totallines;
MyPolyLine SideLine;
std::vector<GeoEdge> SideGeoLine;
std::vector<SimPoint> CultPoints;
std::vector<SimPoint> TreePoints;
osg::ref_ptr<osg::Geode> CreateRoadWithPolyLine(MyPolyLine *line, float width, int Direction, bool hasLight);
void CreateSideLineWithCurrentRoad(osg::ref_ptr<osg::Geode> gd);
osg::ref_ptr<osg::Geode> MergeRoadFromLAndR(osg::ref_ptr<osg::Geode> gLRoad, osg::ref_ptr<osg::Geode> gRRoad, int Direction);
std::vector<std::string> nameParts;
float GetWidthFromName(std::string matname);
int elemNum = 0;
int currShpType = -1;
	
/************************************************************

   MapRand

***********************************************************/

UINT MapRand(UINT nMax)
{
	int nRand = rand();

	float fMap = (float)(nMax)/RAND_MAX;
	float fRetVal = (float)nRand*fMap + 0.5F;

	return (UINT)fRetVal;
};

/************************************************************

   MapPoint

***********************************************************/

IMPLEMENT_DYNAMIC(CMapPoint,CObject)

CMapPoint::CMapPoint()
{
    m_bStatus = 0;
	m_uiIndex = 0;
	m_dbX = 0.0;
	m_dbY = 0.0;
}

CMapPoint::CMapPoint(CMapPoint& pt)
{
   	m_bStatus = pt.m_bStatus;
	m_uiIndex = pt.m_uiIndex;
	m_dbX = pt.m_dbX;
	m_dbY = pt.m_dbY ;
}

CMapPoint::~CMapPoint()
{
}
/*****************************************************************************
  描述:   两点距离
  参数: 　点对象
  返回值  距离

******************************************************************************/
double CMapPoint::Distance(CMapPoint& pt )
{
	return(sqrt((pt.GetX()-m_dbX)*(pt.GetX()-m_dbX)+(pt.GetY()-m_dbY)*(pt.GetY()-m_dbY)));
}

bool  CMapPoint::IsEqual(CMapPoint& pt )
{
	if ( fabs(m_dbX-pt.GetX()) < EP && fabs(m_dbY-pt.GetY()) < EP )
		return TRUE;
	else
		return FALSE;
}

/*****************************************************************************
  描述: 判断指定点是否在线段上
  参数: 　p1 --- 线段起点　p2 --- 线段终点
  返回值  在线段上 返回TRUE 否则返回FALSE

******************************************************************************/

bool   CMapPoint::IsPointInLine(CMapPoint& p1 , CMapPoint& p2 )
{
	double dblMulti;
	// 判断点是否在外接矩形范围内
	if ( m_dbX < min(p1.GetX() ,p2.GetX()) || m_dbX > max(p1.GetX() ,p2.GetX())
	 || m_dbY <  min(p1.GetY() ,p2.GetY()) || m_dbY > max(p1.GetY() ,p2.GetY()))
		return FALSE;
    //计算叉积
	dblMulti = (double)((m_dbX -p1.GetX())*(p2.GetY() -p1.GetY())-(p2.GetX()-p1.GetX())*(m_dbY -p1.GetY()));
    if ( dblMulti == 0 )
		return TRUE;
	else
		return FALSE;
}

/************************************************************

  MapRectangle

***********************************************************/
CMapRectangle::CMapRectangle()
{
	m_dbLeft = 0.0;
	m_dbRight = 0.0;
	m_dbTop = 0.0;
	m_dbBottom = 0.0;
}


CMapRectangle::CMapRectangle(CMapRectangle& MapRectangle )
{
	m_dbLeft = MapRectangle.m_dbLeft;
	m_dbRight = MapRectangle.m_dbRight;
	m_dbTop = MapRectangle.m_dbTop;
	m_dbBottom = MapRectangle.m_dbBottom;
}

CMapRectangle::~CMapRectangle()
{
}

BOOL CMapRectangle::IsPointIn(CMapPoint& point)
{
    if ( min(m_dbLeft,m_dbRight) < point.GetX()  && point.GetX()<max(m_dbLeft,m_dbRight)
		 && min(m_dbTop,m_dbBottom) < point.GetX()  && point.GetY()<max(m_dbTop,m_dbBottom))
	    return TRUE;
	else
		return FALSE;
}

BOOL CMapRectangle::IsInsercet(CMapRectangle& rc)
{
	if (m_dbRight < rc.GetLeft() || m_dbLeft > rc.GetRight()
		|| m_dbTop > rc.GetBottom() || m_dbBottom < rc.GetTop() )
		return FALSE;
	else
		return TRUE;
}

/************************************************************

  MapPoints

***********************************************************/
IMPLEMENT_DYNAMIC(CMapPoints,CObject)

CMapPoints::CMapPoints()
{
}

CMapPoints::CMapPoints(CMapPoints& points)
{
   int i,iCount;
   CMapPoint *pPoint;

   iCount = points.GetCount() - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
	  pPoint = new CMapPoint(*points.Get(i));
	  Add(pPoint);
   }
}

CMapPoints::~CMapPoints()
{
   Clear();
}

long CMapPoints::GetCount()
{
   return m_Points.GetSize();
}

CMapRectangle CMapPoints::GetExtent()
{
	return m_Rectangle;
}


void CMapPoints::SetExtent(CMapRectangle& exent)
{
	m_Rectangle.SetLeft( exent.GetLeft());
	m_Rectangle.SetRight( exent.GetRight());
	m_Rectangle.SetTop(exent.GetTop());
    m_Rectangle.SetBottom(exent.GetBottom());
}

CMapPoint* CMapPoints::Get(long lIndex)
{
	int iCount;
	CMapPoint  *pPt = NULL;

	iCount = m_Points.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return pPt;
    pPt = m_Points.GetAt(lIndex);
	return pPt;
}

void CMapPoints::Add(CMapPoint* pPoint)
{
	if ( pPoint == NULL )
		return;
	m_Points.Add( pPoint );
}

void CMapPoints::Set(long lIndex , CMapPoint* pPoint)
{
	int iCount;
	CMapPoint *pPt;

	if ( pPoint == NULL )
		return;
	iCount = m_Points.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
    pPt = m_Points.GetAt( lIndex );
    m_Points.SetAt(lIndex,pPoint);
    delete pPt;
}

void CMapPoints::Remove(long lIndex)
{
	int iCount;
	CMapPoint *pPt;

	iCount = m_Points.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	pPt = m_Points.GetAt( lIndex );
    m_Points.RemoveAt(lIndex,1);
	delete pPt;
}

void CMapPoints::Insert(long lIndex , CMapPoint* pPoint)
{
	int iCount;

	iCount = m_Points.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	m_Points.InsertAt(lIndex,pPoint);
}

void CMapPoints::Clear()
{
   int i,iCount;
   CMapPoint *pPoint;

   iCount = m_Points.GetSize() - 1;
   for ( i = iCount ; i >= 0   ; i-- )
   {
		pPoint = m_Points.GetAt(i);
		delete pPoint;
   }
   m_Points.RemoveAll();
}

/************************************************************

  MapParts

***********************************************************/
IMPLEMENT_DYNAMIC(CMapParts,CObject)
CMapParts::CMapParts()
{
}

CMapParts::CMapParts(CMapParts& Parts)
{
   int i,iCount;
   CMapPoints *pPoints;

   iCount = Parts.GetCount() - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
		pPoints = new CMapPoints(*Parts.Get(i));
		m_Parts.Add(pPoints);

   }
}

CMapParts::~CMapParts()
{
	Clear();
}

long CMapParts::GetCount()
{
	return m_Parts.GetSize();
}

void CMapParts::Add(CMapPoints* pPoints)
{
	if ( pPoints == NULL )
		return;
	m_Parts.Add(pPoints);

}

void CMapParts::Set(long lIndex, CMapPoints* pPoints)
{
	long lCount;
	CMapPoints* pOldPoints;

	if ( pPoints == NULL )
		return;

	lCount  = m_Parts.GetSize() - 1;
	if ( lIndex < 0 || lIndex > lCount )
		return;
	pOldPoints = m_Parts.GetAt(lIndex);
	m_Parts.SetAt(lIndex , pPoints);
	delete pOldPoints;
}

void CMapParts::Remove(long lIndex)
{
	long lCount;
	CMapPoints *pPoints;

	lCount  = m_Parts.GetSize() - 1;
	if ( lIndex < 0 || lIndex > lCount )
		return;
	pPoints = m_Parts.GetAt(lIndex);
    m_Parts.RemoveAt(lIndex,1);
	delete pPoints;
}

void CMapParts::Insert(long lIndex, CMapPoints* pPoints)
{
	int iCount;

	iCount = m_Parts.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	m_Parts.InsertAt(lIndex,pPoints);
}

CMapPoints* CMapParts::Get(long lIndex)
{

	int iCount;
    CMapPoints  *pPts=NULL;

	iCount = m_Parts.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return pPts;
    pPts = m_Parts.GetAt(lIndex);
	return pPts;
}

void CMapParts::Clear()
{
	long i,lCount;
	CMapPoints* pPoints;

	lCount  = m_Parts.GetSize() - 1;
	for(i = lCount ; i >= 0 ; i-- )
	{
		pPoints = m_Parts.GetAt(i);
		delete pPoints;
    }
	m_Parts.RemoveAll();
}

/************************************************************

  MapLine

***********************************************************/
IMPLEMENT_DYNAMIC(CMapLine,CObject)
CMapLine::CMapLine()
{
}

CMapLine::CMapLine(CMapLine& mapline )
{
   int i,iCount;
   CMapParts *pParts;

   iCount = m_Line.GetSize() - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
		pParts = new CMapParts(*(mapline.GetParts(i)));
		m_Line.Add(pParts);
   }
}

CMapLine::~CMapLine()
{
     Clear();
}

long CMapLine::GetCount()
{
	return m_Line.GetSize();
}

CMapRectangle CMapLine::GetExtent()
{
    return m_Extent;
}

void CMapLine::SetExtent(CMapRectangle& extent)
{
	m_Extent.SetLeft( extent.GetLeft());
	m_Extent.SetRight( extent.GetRight());
	m_Extent.SetTop(extent.GetTop());
    m_Extent.SetBottom(extent.GetBottom());
}

CMapParts* CMapLine::GetParts(long lIndex)
{
	int iCount;
	CMapParts  *pParts = NULL;

	iCount = m_Line.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return pParts;
    pParts = m_Line.GetAt(lIndex);
	return pParts;
}

double CMapLine::GetLength()
{
	return 0;
}
/*****************************************************************************
  描述:   点到线距离
  参数: 　点
  返回值  距离

******************************************************************************/
double CMapLine::Distance(CMapPoint& pt )
{
	int i,j,k;
	double dblMinDist,dblDist;
    CMapParts  *pParts;
	CMapPoints *pPoints;
	CMapPoint  *pPoint1,*pPoint2;

	dblMinDist =  infinity;
	dblDist =  infinity;
	for ( i = 0 ; i < m_Line.GetSize() ; i++ )
	{
		pParts = (CMapParts*)m_Line.GetAt(i);
		for ( j = 0 ; j < pParts->GetCount() ; j++)
		{
			pPoints = (CMapPoints*)pParts;
			for( k = 0 ; k < pPoints->GetCount() - 1 ; k++)
			{
				pPoint1 = (CMapPoint*)pPoints->Get(k);
				pPoint2 = (CMapPoint*)pPoints->Get(k+1);
				//计算点到线段最小距离
				dblDist = ptToSegment(pt,*pPoint1,*pPoint2);
                if ( dblDist <= EP )
					return 0.0;
				else if ( dblDist < dblMinDist )
                    dblMinDist = dblDist;
            }
		}
    }
	return dblMinDist;
}
/*****************************************************************************
  描述:    计算点到线段最小距离
  参数: 　 p1 --- 线段起点　p2 --- 线段终点
  返回值:  在点到线段最小距离

******************************************************************************/
double CMapLine::ptToSegment(CMapPoint& pt,CMapPoint& ptStart,CMapPoint& ptEnd)
{
	double dblDist,dblX,dblY,k;
	CMapPoint pPlumb;

	if ( pt.IsPointInLine(ptStart, ptEnd) )
		return 0.0; //点在直线上
    if ( fabs(ptEnd.GetX() - ptStart.GetX()) <= EP )
	{
		dblX = ptStart.GetX();
		dblY = pt.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.GetY() > min(ptStart.GetY(), ptEnd.GetY()) && pPlumb.GetY()
			< max(ptStart.GetY(), ptEnd.GetY()))
        {
			//垂足在线段范围内
			return fabs(pPlumb.GetX() -pt.GetX());
        }
		else
        {   //判断线段的哪个端点离垂足比较近,然后计算两点间距离
            if ( fabs(pPlumb.GetY() - ptStart.GetY()) < fabs(pPlumb.GetY() - ptEnd.GetY()))
			   dblDist = pt.Distance(ptStart);
			else
			   dblDist = pt.Distance(ptEnd);
        }
		//线段为垂直情况
    } else if ( fabs(ptEnd.GetY() - ptStart.GetY()) <= EP )
    {
		//线段为水平情况
		dblX = pt.GetX();
		dblY = ptStart.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.GetX() > min(ptStart.GetX(), ptEnd.GetX()) && pPlumb.GetX()
			< max(ptStart.GetX(), ptEnd.GetX()))
        {
			//垂足在线段范围内
			return fabs(pPlumb.GetY() -pt.GetY());
        }
		else
        {   //判断线段的哪个端点离垂足比较近,然后计算两点间距离
            if ( fabs(pPlumb.GetX() - ptStart.GetX()) < fabs(pPlumb.GetX() - ptEnd.GetX()))
			   dblDist = pt.Distance(ptStart);
			else
			   dblDist = pt.Distance(ptEnd);
        }
    }
	else
    {
		//线段倾斜状态
		k = (double)((ptEnd.GetY() - ptStart.GetY() ) /(ptEnd.GetX() - ptStart.GetX()));
		// 计算垂足坐标
		dblX = (k*k*ptStart.GetX() + k*(pt.GetY()-ptStart.GetY())+pt.GetX())/(k*k+1);
		dblY = k*(dblX-ptStart.GetX()) + ptStart.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.IsPointInLine(ptStart,ptEnd) )
        {
			//垂足在线段上
		    dblDist = pt.Distance(pPlumb);
        }
		else
        {
			//判断线段的哪个端点离垂足比较近,然后计算两点间距离
			if ( pPlumb.Distance(ptStart) < pPlumb.Distance(ptEnd))
               dblDist = pt.Distance(ptStart);
			else
			   dblDist = pt.Distance(ptEnd);
		}
    }
	return dblDist;
}

void CMapLine::Add(CMapParts* pParts)
{

    if ( pParts == NULL )
		return;
	m_Line.Add( pParts );
}

void CMapLine::Set(long lIndex , CMapParts* pParts)
{
	int iCount;
	CMapParts *pOldParts;

	if ( pParts == NULL )
		return;

	iCount = m_Line.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
    pOldParts = m_Line.GetAt( lIndex );
    m_Line.SetAt(lIndex,pParts);
    delete pOldParts;
}

void CMapLine::Remove(long lIndex)
{
	int iCount;
	CMapParts *pParts;

	iCount = m_Line.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	pParts = m_Line.GetAt( lIndex );
    m_Line.RemoveAt(lIndex,1);
	delete pParts;
}

void CMapLine::Insert(long lIndex , CMapParts* pParts)
{
	int iCount;

	iCount = m_Line.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	m_Line.InsertAt(lIndex,pParts);
}

void CMapLine::Clear()
{
   int i,iCount;
   CMapParts *pParts;

   iCount = m_Line.GetSize() - 1;
   for ( i = iCount ; i >= 0   ; i-- )
   {
		pParts = m_Line.GetAt(i);
		delete pParts;
   }
   m_Line.RemoveAll();
}

/************************************************************

  MapPolygon

***********************************************************/
IMPLEMENT_DYNAMIC(CMapPolygon,CObject)
CMapPolygon::CMapPolygon()
{
}

CMapPolygon::CMapPolygon(CMapPolygon& mappolygon )
{
   int i,iCount;
   CMapParts *pParts;

   iCount = m_Polygon.GetSize() - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
		pParts = new CMapParts(*mappolygon.GetParts(i));
		m_Polygon.Add(pParts);
   }
}

CMapPolygon::~CMapPolygon()
{
     Clear();
}

long CMapPolygon::GetCount()
{
	return m_Polygon.GetSize();
}

CMapRectangle CMapPolygon::GetExtent()
{
    return m_Extent;
}

void CMapPolygon::SetExtent(CMapRectangle& extent)
{
	m_Extent.SetLeft( extent.GetLeft());
	m_Extent.SetRight( extent.GetRight());
	m_Extent.SetTop(extent.GetTop());
    m_Extent.SetBottom(extent.GetBottom());
}

CMapParts* CMapPolygon::GetParts(long lIndex)
{
	int iCount;
    CMapParts  *pParts;

	iCount = m_Polygon.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return NULL;
    pParts = m_Polygon.GetAt(lIndex);
	return pParts;
}

double CMapPolygon::GetArea()
{
	return 0;
}

void CMapPolygon::Add(CMapParts* pParts)
{
    if ( pParts == NULL )
		return;
	m_Polygon.Add( pParts );
}

void CMapPolygon::Set(long lIndex , CMapParts* pParts)
{
	int iCount;
	CMapParts *pOldParts;

	if ( pParts == NULL )
		return;

	iCount = m_Polygon.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
    pOldParts = m_Polygon.GetAt( lIndex );
    m_Polygon.SetAt(lIndex,pParts);
    delete pOldParts;
}

void CMapPolygon::Remove(long lIndex)
{
	int iCount;
	CMapParts *pParts;

	iCount = m_Polygon.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	pParts = m_Polygon.GetAt( lIndex );
    m_Polygon.RemoveAt(lIndex,1);
	delete pParts;
}

void CMapPolygon::Insert(long lIndex , CMapParts* pParts)
{
	int iCount;

	iCount = m_Polygon.GetSize()-1;
	if ( lIndex < 0 || lIndex > iCount )
		return ;
	m_Polygon.InsertAt(lIndex,pParts);
}

void CMapPolygon::Clear()
{

   int i,iCount;
   CMapParts *pParts;

   iCount = m_Polygon.GetSize() - 1;
   for ( i = iCount ; i >= 0   ; i-- )
   {
		pParts = m_Polygon.GetAt(i);
		delete pParts;
   }
   m_Polygon.RemoveAll();
}

/*************************************************
  描述:         判断点是否在多边形内
  输入:         点对象
  输出：        在返回TRUE 不在返回FALSE
*************************************************/
BOOL CMapPolygon::IsPointIn(CMapPoint& pt )
{
    int i,j,k,iNumber;
	double dblTemp;
	CMapParts*  pParts;
	CMapPoints* pPoints;
	CMapPoint  *ptFirst,*ptSecond;
	CMapPoint	ptInfint;

	iNumber = 0;
	//不在矩形范围内直接返回
	if ( !m_Extent.IsPointIn(pt) )
		  return FALSE;
	//做一条通过pt点的射线
	dblTemp = pt.GetY();
	ptInfint.SetY( dblTemp );
	dblTemp = infinity;
	ptInfint.SetX(dblTemp);
	//获得复合多边形的每一个多边形
	for ( i = 0 ; i < m_Polygon.GetSize() ; i++ )
    {
		pParts = (CMapParts*)m_Polygon.GetAt(i);
		//获得一个多边形的点集合
		for ( j = 0 ; j < pParts->GetCount() ; j++ )
		{
			pPoints = (CMapPoints*)pParts->Get(j);
			//获得每个点集合的顶点坐标
			for ( k = 0 ; k < pPoints->GetCount() - 1 ; k++ )
			{
				ptFirst  = (CMapPoint*)pPoints->Get(k);
				ptSecond = (CMapPoint*)pPoints->Get(k+1);
				if (pt.IsPointInLine(*ptFirst,*ptSecond) )
					return TRUE; //该点在多边形边上
                if ( ptSecond->GetY() == ptFirst->GetY() )
					continue;  //略过水平边

				if ( ptFirst->IsPointInLine(  pt , ptInfint ))
				{
					//边与射线相交于边第一个顶点
					if ( ptFirst->GetY() > ptSecond->GetY() )
						iNumber++;
				}
				else if( ptSecond->IsPointInLine( pt , ptInfint ))
				{
					//边与射线相交于边的第二个顶点
					if ( ptFirst->GetY() > ptSecond->GetY() )
					    iNumber++;

				} //判断是否有公共交点
				else if ( isIntersect(pt,ptInfint,*ptFirst,*ptSecond) )
						iNumber++ ;
			}
		}

	}
	if ( iNumber % 2 )
		return TRUE;
    else
	    return FALSE;
}

/*****************************************************************************
  描述:  判断两条直线是否相交 且交点不是顶点
  参数: p1 --- 线段起点 p2 ---线段终点 　p3 --- 线段起点　p4 --- 线段终点
  返回值  在直线上 返回TRUE 否则返回FALSE
******************************************************************************/

BOOL CMapPolygon::isIntersect(CMapPoint& p1 , CMapPoint& p2 , CMapPoint& p3 , CMapPoint& p4 )
{
	double dblMulti,dblTmp1,dblTmp2;
	//顶点相交
	if ( p1.IsEqual(p3) || p1.IsEqual(p4) || p2.IsEqual(p3) || p2.IsEqual(p4) )
		return FALSE;
    //判断两条线段外接矩形是否相交
	if ( min(p1.GetX(),p2.GetX()) > max(p3.GetX(),p4.GetX()) || max(p1.GetX(),p2.GetX())
		< min(p3.GetX(),p4.GetX()) || min(p1.GetY(),p2.GetY()) > max(p3.GetY(),p4.GetY())
		|| max(p1.GetY(),p2.GetY()) < min(p3.GetY(),p4.GetY()))
		return FALSE;
    //计算叉积
	dblTmp1 = (double)((p1.GetX() - p3.GetX())*(p4.GetY()-p3.GetY()) - (p4.GetX()-p3.GetX())*(p1.GetY() - p3.GetY()));
	dblTmp2 = (double)((p4.GetX() -p3.GetX())*(p2.GetY() - p3.GetY()) - (p2.GetX()-p3.GetX())*(p4.GetY()-p3.GetY()));        ;
	dblMulti = dblTmp1 * dblTmp2;

	if ( dblMulti >= 0 )
		return TRUE;
	else
		return false;
}

/************************************************************

  MapField

***********************************************************/
CMapField::CMapField()
{
	m_csFieldName = _T("");
	m_csValue  = _T("");
	m_lFieldType = fdInvaild;
	::VariantInit(&m_varValue);
}

CMapField::CMapField(CMapField& field)
{
	m_csFieldName = field.GetName();
	m_csValue  = field.GetValueAsString();
	m_lFieldType = field.GetType();
	m_varValue = field.GetValue();
}

CMapField::~CMapField()
{
}

CString CMapField::GetName()
{
	return m_csFieldName;
}

void CMapField::SetName(LPCTSTR lpszName)
{
	m_csFieldName = lpszName;
}

long CMapField::GetType()
{
	return m_lFieldType;
}

void CMapField::SetType(long lType)
{
	m_lFieldType = lType;
}

CString CMapField::GetValueAsString()
{
	CString csValue = _T("");
	switch( m_lFieldType )
    {
		case fdInteger:
			csValue.Format("%d",m_varValue.lVal);
			break;
        case fdDouble:
			csValue.Format("%f",m_varValue.dblVal);
			break;
        case fdString:
			return m_csValue;
		    break;
		case fdInvaild:
			break;
		default:
			break;
    }
	return csValue;
}

void CMapField::SetValueAsString(LPCTSTR lpstr)
{
	m_csValue = lpstr;
}

VARIANT CMapField::GetValue()
{
	return m_varValue;
}

void CMapField::SetValue(const VARIANT& var)
{
    switch( m_lFieldType )
    {
		case fdInteger:
		    m_varValue.bVal = var.bVal;
			m_varValue.lVal = var.lVal;
			break;
        case fdDouble:
			m_varValue.bVal = var.bVal;
			m_varValue.dblVal = var.dblVal;
			break;
       	case fdInvaild:
			break;
		default:
			m_varValue = var;
			break;
    }

}

/************************************************************

  MapFields

***********************************************************/

CMapFields::CMapFields()
{
}

CMapFields::CMapFields(CMapFields& fields )
{
   int i,iCount;

   CMapField *pField;

   iCount = m_Fields.GetSize() - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
		pField = new CMapField(*(fields.GetField(i)));
		m_Fields.Add(pField);
   }
}

CMapFields::~CMapFields()
{
	Clear();
}


short CMapFields::GetCout()
{

	return m_Fields.GetSize();
}

void CMapFields::Add(CMapField* pField)
{
	if ( pField == NULL )
		return;
	m_Fields.Add( pField );
}


void CMapFields::Remove(short sIndex)
{
	int iCount;
	CMapField *pField;

	iCount = m_Fields.GetSize()-1;
	if ( sIndex < 0 || sIndex > iCount )
		return ;
	pField = m_Fields.GetAt( sIndex );
    m_Fields.RemoveAt(sIndex,1);
	delete pField;
}

void CMapFields::Insert(short sIndex, CMapField* pField)
{

	int iCount;

	iCount = m_Fields.GetSize()-1;
	if ( sIndex < 0 || sIndex > iCount )
		return ;
	m_Fields.InsertAt(sIndex,pField);
}

CMapField* CMapFields::GetField(short sIndex)
{
	int iCount;
	CMapField  *pField = NULL;

	iCount = m_Fields.GetSize()-1;
	if ( sIndex < 0 || sIndex > iCount )
		return pField;
    pField = m_Fields.GetAt(sIndex);
	return pField;
}

void CMapFields::Clear()
{
	long i,lCount;
	CMapField* pField;

	lCount  = m_Fields.GetSize() - 1;
	for(i = lCount ; i >= 0 ; i-- )
	{
		pField = m_Fields.GetAt(i);
		delete pField;
    }
	m_Fields.RemoveAll();
}




/************************************************************

  MapTableDesc

***********************************************************/

CMapTableDesc::CMapTableDesc()
{


}

CMapTableDesc::CMapTableDesc(CMapTableDesc& tblDesc)
{

   int i,iCount;
   FIELD_ELEMENT *pField,*pSource;

   iCount = tblDesc.GetFieldCount()  - 1;
   for ( i = 0 ; i <= iCount ; i++ )
   {
		pField = new  FIELD_ELEMENT;
		pSource = tblDesc.GetDesc(i);
        strcpy(pField->szFieldName,pSource->szFieldName);
		pField->cFieldType = pSource->cFieldType;
		pField->ucFieldDecimal = pSource->ucFieldDecimal;
		pField->ucFieldLength = pSource->ucFieldLength;
		pField->ulOffset = pSource->ulOffset;
		pField->cProductionIndex = pSource->cProductionIndex;
		pField->dbaseiv_id = pSource->dbaseiv_id;
		strcpy(pField->reserved1 , pField->reserved1);
		strcpy(pField->reserved2 , pField->reserved2);
		m_fieldsDesc.Add(pField);
   }
}
CMapTableDesc::~CMapTableDesc()
{
	Clear();

}

short CMapTableDesc::GetFieldCount()
{
	return m_fieldsDesc.GetSize();

}

void CMapTableDesc::SetFieldCount(short sCount)
{


}

CString CMapTableDesc::GetFieldName(short sIndex)
{
	int iCount;
	FIELD_ELEMENT *pField;

	iCount = m_fieldsDesc.GetSize()  - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return _T("");
    pField = m_fieldsDesc.GetAt(sIndex);
	CString csName(pField->szFieldName);
	return csName;
}

void CMapTableDesc::SetFieldName(short sIndex, LPCTSTR lpszNewValue)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return ;
    pField = m_fieldsDesc.GetAt(sIndex);
	strcpy(pField->szFieldName , lpszNewValue);

}

long CMapTableDesc::GetFieldType(short sIndex)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize()  - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return INVALID_VALUE;
    pField = m_fieldsDesc.GetAt(sIndex);
	return pField->cFieldType;

}

void CMapTableDesc::SetFieldType(short sIndex, long nNewValue)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return;
    pField = m_fieldsDesc.GetAt(sIndex);
	pField->cFieldType  = (unsigned char )nNewValue;


}

short CMapTableDesc::GetFieldPrecision(short sIndex)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return INVALID_VALUE;
    pField = m_fieldsDesc.GetAt(sIndex);
	return pField->ucFieldDecimal;

}

void CMapTableDesc::SetFieldPrecision(short sIndex, short nNewValue)
{
    int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return;
    pField = m_fieldsDesc.GetAt(sIndex);
	pField->ucFieldDecimal  = (unsigned char)nNewValue;

}

short CMapTableDesc::GetFieldLength(short sIndex)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return INVALID_VALUE;
    pField = m_fieldsDesc.GetAt(sIndex);
	return pField->ucFieldLength;
}

void CMapTableDesc::SetFieldLength(short sIndex, short nNewValue)
{

    int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return;
    pField = m_fieldsDesc.GetAt(sIndex);
	pField->ucFieldLength  = (unsigned char)nNewValue;
}

/*short CMapTableDesc::GetFieldScale(short sIndex)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return INVALID_VALUE;
    pField = m_fieldsDesc.GetAt(sIndex);
	return pField->sfdScale;

}

void CMapTableDesc::SetFieldScale(short sIndex, short nNewValue)
{
	int iCount;
	FIELD_ELEMENT* pField;

	iCount = m_fieldsDesc.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount )
		return;
    pField = m_fieldsDesc.GetAt(sIndex);
	pField->sfdScale  = nNewValue;

}*/

FIELD_ELEMENT* CMapTableDesc::GetDesc(short sIndex)
{
	int iCount;
	FIELD_ELEMENT*  pDesc = NULL;

	iCount = m_fieldsDesc.GetSize()-1;
	if ( sIndex < 0 || sIndex > iCount )
		return pDesc;
    pDesc = m_fieldsDesc.GetAt(sIndex);
	return pDesc;
}
void CMapTableDesc::Add(FIELD_ELEMENT* pfdDesc)
{

    if ( pfdDesc == NULL )
		return;
	m_fieldsDesc.Add( pfdDesc);
}

void CMapTableDesc::Remove(short sIndex)
{
	int iCount;
	FIELD_ELEMENT* pfdDesc;

	iCount = m_fieldsDesc.GetSize()-1;
	if ( sIndex < 0 || sIndex > iCount )
		return ;
	pfdDesc = m_fieldsDesc.GetAt( sIndex );
    m_fieldsDesc.RemoveAt(sIndex,1);
	delete pfdDesc;
}

void CMapTableDesc::Insert(short sindex,FIELD_ELEMENT* pfdDesc)
{
	int iCount;

	iCount = m_fieldsDesc.GetSize()-1;
	if ( sindex < 0 || sindex > iCount )
		return ;
	m_fieldsDesc.InsertAt(sindex,pfdDesc);

}

void CMapTableDesc::Clear()
{
   int i,iCount;
   FIELD_ELEMENT *pfdDesc;

   iCount = m_fieldsDesc.GetSize() - 1;
   for ( i = iCount ; i >= 0   ; i-- )
   {
		pfdDesc = m_fieldsDesc.GetAt(i);
		delete pfdDesc;
   }
   m_fieldsDesc.RemoveAll();
}

/************************************************************

  MapRecordSet

***********************************************************/

CMapRecordSet::CMapRecordSet()
{

	m_CacheSize = 50;
	bBOF = FALSE;
	bEOF = FALSE;
	m_bDbfOpen = FALSE;
	iCursorPos = -1;
}

CMapRecordSet::~CMapRecordSet()
{

	Clear();
}

long CMapRecordSet::GetRecordCount()
{

	return m_Header.no_recs;
}

CMapFields* CMapRecordSet::GetFields(long lIndex)
{
	long lCount;
    CMapFields *fields = NULL;
	CMapFields *pFields;
	lCount = m_Fields.GetSize() - 1;
    if ( lIndex < 0 || lIndex > lCount )
		return fields;
	pFields = m_Fields.GetAt(lIndex);
	return pFields;
}

CMapTableDesc* CMapRecordSet::GetTableDesc()
{

	return &m_TableDesc;
}

BOOL CMapRecordSet::GetBOF()
{
	return bBOF;

}

BOOL CMapRecordSet::GetEOF()
{
	return bEOF;

}

int  CMapRecordSet::GetCacheSize()
{
	return m_CacheSize;

}

BOOL CMapRecordSet::SetCacheSize(int& CacheSize)
{
	if ( CacheSize < 0 || CacheSize > MAX_CACH_SIZE )
		return FALSE;
	return TRUE;
}

/*************************************************
  描述:         打开DBF文件
  输入:         文件名
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/
BOOL CMapRecordSet::openDBF(CString& csFileName)
{
	unsigned int  iTemp;
	unsigned long ulReocrdCount;
	unsigned short ulLength,ulRecLength;
	short i,sFieldCount;
	char*  pszBuffer;
	FIELD_ELEMENT *pField,*pOldField;
    CFileException fe;

	//打开主文件
	if ( !fDbf.Open(csFileName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	m_bDbfOpen = TRUE;
	fDbf.Seek(0L, CFile::begin);
	//读入文件头
	if ( fDbf.Read(&m_Header,sizeof(m_Header)) != sizeof(m_Header))
    	return FALSE;
	ulReocrdCount = m_Header.no_recs;
	ulLength = m_Header.head_len;
    ulRecLength = m_Header.rec_len;

	//计算字段个数
	sFieldCount = (ulLength - sizeof(DBF_HEADER)-1)/sizeof(FIELD_ELEMENT);
	iTemp = sFieldCount * sizeof(FIELD_ELEMENT) + 1;
	pszBuffer = new char[iTemp];
	if ( pszBuffer == NULL )
		return FALSE;
    //读入字段描述部分数据(表结构)
	if ( fDbf.Read(pszBuffer,iTemp) != iTemp)
    {
    	delete []pszBuffer;
		return FALSE;
	}
	for ( i = 0 ; i < sFieldCount ; i++ )
	{
		pField = new FIELD_ELEMENT;
		if ( pField == NULL )
		{
			delete []pszBuffer;
		    m_TableDesc.Clear();
			return FALSE;
        }

		memcpy(pField,pszBuffer+i*sizeof(FIELD_ELEMENT),sizeof(FIELD_ELEMENT));
		if ( i == 0 )
        	pField->ulOffset = 0;
        else
			pField->ulOffset = pOldField->ulOffset + pOldField->ucFieldLength;

		//判断字段类型
		if ( pField->cFieldType != 'N' && pField->cFieldType != 'F' )
        {
			pField->ucFieldLength += pField->ucFieldDecimal*256;
            pField->ucFieldDecimal = 0;
		}
		pOldField = pField;
       	m_TableDesc.Add(pField);

    }
    //读入记录到记录集缓冲区
    ReadRecord(0);

	delete []pszBuffer;
	return true;
}

/*************************************************
  描述:         读入数据记录
  输入:         记录的索引(从0开始)
  输出：        无
*************************************************/

void CMapRecordSet::ReadRecord(unsigned long lRecordID)
{
	int	   j,iRecordOffset;
	char   *pszBuffer;
	char   szBuff[255];
	double dbValue;
	FIELD_ELEMENT *pField;
	CMapField *pMapField;
	CMapFields *pMapFields;

	VARIANT varValue;

	if( lRecordID < 0 || lRecordID >= m_Header.no_recs)
        return; //无效的索引值

    if ( iCursorPos != lRecordID) //要读取的记录未在缓存中
    {
		//计算记录相对文件头的偏移量
		iRecordOffset = lRecordID*m_Header.rec_len + m_Header.head_len;
		pszBuffer = new char[ m_Header.rec_len];
        fDbf.Seek(iRecordOffset , CFile::begin);
		if ( fDbf.Read(pszBuffer,m_Header.rec_len) != m_Header.rec_len)
		{
    		delete []pszBuffer;
			return ;
		}

	    Clear();
		pMapFields = new CMapFields();
		for ( j = 0 ; j < m_TableDesc.GetFieldCount() ; j++)
		{
			pMapField = new CMapField;
			pField = m_TableDesc.GetDesc(j);
			pMapField->SetName( pField->szFieldName );
			pMapField->SetType(pField->cFieldType);
			memset(szBuff , 0 , 255);
			//略过该记录是否删除标记字节pszBuffer+1
			strncpy(szBuff, pszBuffer+1+pField->ulOffset , pField->ucFieldLength);
			if ( pField->cFieldType == 'N' || pField->cFieldType == 'F' )
			{
				::VariantInit(&varValue);
				dbValue = atof(szBuff );
				if ( pField->ucFieldDecimal == 0 )
				{
					varValue.bVal = VT_I4;
					varValue.lVal = (int)dbValue;
					pMapField->SetType(fdInteger);
				}
				else
				{
					varValue.bVal = VT_R8;
					varValue.dblVal = dbValue;
					pMapField->SetType(fdDouble);
		   		}
				pMapField->SetValue(varValue);
			}
			else if ( pField->cFieldType == 'C' )
			{
				pMapField->SetValueAsString(szBuff);
				pMapField->SetType(fdString);
			}
			else
				pMapField->SetType(fdInvaild);
			pMapFields->Add(pMapField);

		}
		Clear(); //清空缓冲区加入新记录
		m_Fields.Add(pMapFields);
		delete []pszBuffer;
		iCursorPos = lRecordID;
	}
}

/*************************************************
  描述:         移动到记录集头部
  输入:         无
  输出：        无
*************************************************/
void CMapRecordSet::MoveFirst()
{
	 if ( !m_bDbfOpen ) //数据库文件未打开
		 return;
	 bBOF = TRUE;
	 bEOF = FALSE;
	 iCursorPos = -1;
	 ReadRecord(0);
}

/*************************************************
  描述:         移动到记录集尾部
  输入:         无
  输出：        无
*************************************************/
void CMapRecordSet::MoveLast()
{
     if ( !m_bDbfOpen )
		 return;
	 bEOF = TRUE;
	 ReadRecord(m_Header.no_recs - 1);
}

/*************************************************
  描述:         移动到下一条记录
  输入:         无
  输出：        无
*************************************************/

void CMapRecordSet::MoveNext()
{

	 if ( !m_bDbfOpen )
		 return;

	 if ( m_Header.no_recs == 1 )
     {
		 bBOF = TRUE;
		 bEOF = TRUE;
		 return;
     }
	 if ( iCursorPos < m_Header.no_recs-1)
     {
		 bBOF = FALSE;
		 ReadRecord(iCursorPos + 1);
     }
	 else
        bEOF = TRUE;
}

/*************************************************
  描述:         移动到上一条记录
  输入:         无
  输出：        无
*************************************************/
void CMapRecordSet::MovePrev()
{
	if ( !m_bDbfOpen ) //数据库文件未打开
		 return;

	if ( m_Header.no_recs == 1 )
    {
		 bBOF = TRUE;
		 bEOF = TRUE;
		 return;

    }
	if ( iCursorPos > 0 )
    {
		 bEOF = FALSE;
		 ReadRecord(iCursorPos - 1);

	}
	else
        bBOF = TRUE;
}

/*************************************************
  描述:         移动iNumRecords条记录
  输入:         移动的记录数、移动相对位置
  输出：        无
*************************************************/

BOOL CMapRecordSet::Move(int iNumRecords , RECORDSTART Start )
{
	int iPos;
	/*if ( bEOF && iNumRecords > 0 ) //已经到记录集末尾
		return FALSE;
	if ( bBOF && iNumRecords < 0 ) //已经到记录集头
		return FALSE;
	if ( iNumRecords == 0 )
		return TRUE;*/

	switch ( Start )
    {
		case BookmarkCurrent:
			iPos = iCursorPos;
			break;
        case BookmarkLast:
			iPos = m_Header.no_recs - 1;
			break;
        default: // BookmarkFirst
			iPos = 0;
			break;
	}

	if ( iNumRecords > 0 ) //向后移动
    {
		if ( m_Header.no_recs <= (unsigned long)(iPos + iNumRecords))
			return FALSE;
		else
        {
			ReadRecord(iPos + iNumRecords);
			return TRUE;
		}
	}
	else
    {
	    if (  (iPos + iNumRecords) < 0 )
			return FALSE;
		else
        {
			ReadRecord(iPos + iNumRecords);
			return TRUE;
		}
    }
}

/*************************************************
  描述:         清空记录缓冲区
  输入:         无
  输出：        无
*************************************************/
void CMapRecordSet::Clear()
{
	int i;
	CMapFields *pMapFields;

	for ( i = m_Fields.GetSize() - 1 ; i >=0  ; i-- )
	{
		pMapFields = m_Fields.GetAt(i);
		delete pMapFields;
	}
	m_Fields.RemoveAll();
}



/************************************************************

  MapRender

***********************************************************/
CMapRender::CMapRender()
{
	m_SimpleRender.FillColor = RGB(MapRand(255),MapRand(255),MapRand(255));
	m_SimpleRender.OutlineColor = RGB(0,0,0);
	m_SimpleRender.iIndex = 0;
	m_RenderType = SIMPLE_RENDER;
	m_FieldIndex = -1;
}


CMapRender::~CMapRender()
{
	Clear();
}

void CMapRender::Add(RENDERINFO& rInfo )
{
	int i;
	RENDERINFO* pInfo;

	for ( i = m_Render.GetSize() - 1 ; i >= 0 ; i--)
    {

		pInfo = (RENDERINFO*)m_Render.GetAt(i);
	    if ( pInfo->csValue == rInfo.csValue )
			return;
	}
	pInfo = new RENDERINFO;
    if ( pInfo == NULL )
		return;
	pInfo->csValue = rInfo.csValue;
    pInfo->clr =  rInfo.clr;
	m_Render.Add(pInfo);
}

void  CMapRender::RemoveByValue(CString& csValue)
{
	int i;
	RENDERINFO* pInfo;

	for ( i = m_Render.GetSize() - 1 ; i >= 0 ; i--)
    {
		pInfo = (RENDERINFO*)m_Render.GetAt(i);
	    if ( pInfo->csValue == csValue )
        {
			delete pInfo;
			m_Render.RemoveAt(i);
			break;
        }
    }
}

void  CMapRender::RemoveByIndex(int iIndex )
{
	RENDERINFO* pInfo;

	if ( iIndex < 0 || iIndex >= m_Render.GetSize())
		return ;

	pInfo = (RENDERINFO*)m_Render.GetAt(iIndex);
	delete pInfo;
	m_Render.RemoveAt(iIndex);
}

int  CMapRender::GetCount()
{
	return m_Render.GetSize();
}


RENDERINFO* CMapRender::GetByValue(CString& csValue)
{
	int i;
	RENDERINFO* pInfo;

	for ( i = m_Render.GetSize() - 1 ; i >= 0 ; i--)
    {

		pInfo = (RENDERINFO*)m_Render.GetAt(i);
	    if ( pInfo->csValue == csValue )
			return pInfo;
    }
	return NULL;
}

RENDERINFO* CMapRender::GetByIndex(int iIndex)
{
	RENDERINFO* pInfo;

	if ( iIndex < 0 || iIndex >= m_Render.GetSize())
		return NULL;

	pInfo = (RENDERINFO*)m_Render.GetAt(iIndex);
	return pInfo;
}

void CMapRender::Clear()
{
	int i;
	RENDERINFO* pInfo;

	for ( i = m_Render.GetSize() - 1 ; i >= 0 ; i--)
    {

		pInfo = (RENDERINFO*)m_Render.GetAt(i);
		delete pInfo;
    }
	m_Render.RemoveAll();
}

void CMapRender::SetSimpleRender( SIMPLERENDER& simpleRender )
{
	m_SimpleRender.FillColor = simpleRender.FillColor;
	m_SimpleRender.OutlineColor = simpleRender.OutlineColor;
	m_SimpleRender.iIndex = simpleRender.iIndex;
}

void CMapRender::GetSimpleRender( SIMPLERENDER& simpleRender )
{
	simpleRender.FillColor  =	m_SimpleRender.FillColor;
	simpleRender.OutlineColor = m_SimpleRender.OutlineColor;
	simpleRender.iIndex = m_SimpleRender.iIndex ;
}

void CMapRender::SetFieldIndex(int iIndex)
{
	m_FieldIndex = iIndex;
}

int CMapRender::GetFieldIndex()
{
	return m_FieldIndex ;
}

void CMapRender::Clone(CMapRender *pRender)
{
	int i;
	SIMPLERENDER simple;
	RENDERINFO   *pInfo;
	RENDERINFO   InfoNew;
	if ( pRender == NULL )
		return;
	pRender->SetRenderType(GetRenderType());
    GetSimpleRender(simple);
	pRender->SetSimpleRender(simple);
	pRender->SetFieldIndex(GetFieldIndex());

    for ( i = 0 ; i < m_Render.GetSize() ; i++ )
    {
		pInfo = (RENDERINFO*)m_Render.GetAt(i);
		InfoNew.csValue = pInfo->csValue;
		InfoNew.clr  = pInfo->clr;
        pRender->Add(InfoNew);

    }
}

/************************************************************

  shpFile

***********************************************************/

CShpFile::CShpFile()
{
	int i;
	bShpOpen = FALSE;
	bShxOpen = FALSE;
	m_shpType = NULLSHP;

	/* -------------------------------------------------------------------- */
	/*	Establish the byte order on this machine.			*/
	/* -------------------------------------------------------------------- */
	i = 1;
	if( *((uchar *) &i) == 1 )
		m_bBigEndian = FALSE;
	else
		m_bBigEndian = TRUE;
	isLake = false;
	islandInOcean = false;
	isAirport = false;
	m_elevOfLake = 0.0;
	m_elevOfAirport = 0.0;
}

CShpFile::~CShpFile()
{

	 CMapPoint*    pPt;
	 CMapPoints*   pPts;
	 CMapLine*     pLine;
	 CMapPolygon*  pPolygon;

	 while(!m_ObList.IsEmpty())
	 {
		switch ( m_shpType )
        {
			case POINT:
				pPt =(CMapPoint*)m_ObList.RemoveTail();
				delete pPt;
				break;
            case MULTIPOINT:
				pPts =(CMapPoints*)m_ObList.RemoveTail();
				delete pPts;
				break;

			case POLYLINE:
                pLine =(CMapLine*)m_ObList.RemoveTail();
				delete pLine;
				break;
			case POLYGON:
                pPolygon =(CMapPolygon*)m_ObList.RemoveTail();
				delete pPolygon;
				break;
            default:
				m_ObList.RemoveTail();
				break;

		}
     }
	 if ( bShpOpen )
		fShp.Close();
	 if (bShxOpen)
		fShx.Close();
	 if ( pRecordSet )
		 delete pRecordSet;

}

/*************************************************
  描述:         读入Shp、shx文件 格式 参阅ESRI技术白皮书
  输入:         文件名(全路径)
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/


BOOL CShpFile::ReadShp(CString& csFileName)
{
	int	  iTemp;
	CString csShxName;
	CFileException fe;
    SHPHEADER varHeader;

	//打开主文件
	if ( !fShp.Open(csFileName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	bShpOpen = TRUE;

	//打开索引文件
	csShxName = csFileName.Left(csFileName.GetLength() - 3);

	csShxName = csShxName + "shx";
    if ( !fShx.Open(csShxName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	bShxOpen = TRUE;

	TRY
    {
	    //读主文件头 长100字节
		if ( fShp.Read(&varHeader , sizeof(SHPHEADER))!= sizeof(SHPHEADER))
			return FILE_READERR;

        iTemp = varHeader.iFileCode;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&iTemp);
		if ( iTemp != 9994 ) //是否是shp文件
			return FILE_CODEERR;
        if ( varHeader.iVersion != FILE_VERSION ) //文件版本是否正确
			return FILE_VERSIONERR;

		//shp 类型
		m_shpType = varHeader.iShpType;
		m_shpFileLength = varHeader.iFileLength;
	    if ( !m_bBigEndian )
			SwapWord(sizeof(int),&m_shpFileLength);

		//保存数据最大矩形范围
		m_Extent.SetLeft(varHeader.dbXMin);
		m_Extent.SetRight(varHeader.dbXMax);
		m_Extent.SetTop(varHeader.dbYMin);
		m_Extent.SetBottom(varHeader.dbYMax);

		//读索引文件头 长100字节
		if ( fShx.Read(&varHeader , sizeof(SHPHEADER))!= sizeof(SHPHEADER))
			return FILE_READERR;
		iTemp = varHeader.iFileCode;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&iTemp);
		if ( iTemp != 9994 ) //是否是shx文件
			return FILE_CODEERR;
        if ( varHeader.iVersion != FILE_VERSION ) //文件版本是否正确
			return FILE_VERSIONERR;
		m_shxFileLength = varHeader.iFileLength;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&m_shxFileLength);
		//通过索引文件计算主文件记录个数 文件长度数值以16位计
		m_iRecordCount = ((m_shxFileLength - 50 )*2)/sizeof(SHXRECORD);

		if ( !ReadRecord() )
			return FILE_READERR;
         if ( !ReadDBF(csFileName))
			 return FALSE;
	}
	CATCH(CFileException ,eload)
    {
		fShp.Abort();
		return FALSE;
    }
	END_CATCH


	return TRUE;

}
/*************************************************
  描述:         读入DBF文件格式
  输入:         文件名(全路径)
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/
BOOL CShpFile::ReadDBF(CString& csFileName)
{
	CString csDbfName;
	BOOL bResult;
	//创建记录集对象
	pRecordSet = new CMapRecordSet;
	ASSERT ( pRecordSet != NULL );
	csDbfName = csFileName.Left(csFileName.GetLength() - 3);
	csDbfName = csDbfName + "dbf";
	//打开DBF文件
	bResult =  pRecordSet->openDBF(csDbfName);
    if ( !bResult )
		delete pRecordSet;
	return bResult;
}

/*************************************************
  描述:         获得当前地图文件最大矩形范围
  输入:         无
  输出：        地图矩形范围
*************************************************/

CMapRectangle CShpFile::GetExtent()
{
	return m_Extent;
}

/*************************************************
  描述:         设置地图文件最大矩形范围
  输入:         地图矩形范围
  输出：        无
*************************************************/
void CShpFile::SetExtent(CMapRectangle& extent )
{
	m_Extent.SetLeft(extent.GetLeft());
    m_Extent.SetRight(extent.GetRight());
	m_Extent.SetTop(extent.GetTop());
    m_Extent.SetBottom(extent.GetBottom());
}

/*************************************************
  描述:         设置Shp对象类型
  输入:         shp 类型
  输出：        无
*************************************************/
void CShpFile::SetShpType(int& iShpType)
{
	m_shpType = iShpType;
}

/*************************************************
  描述:         读入Shp对象坐标
  输入:         无
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/

BOOL CShpFile::ReadRecord()
{
	int i,j,k,m;
	int iDataLen,iLength,iIndex;
	int *pIParts;
	double dbTmp;
	char *pszBuffer,*lpVer;
	BOOL bEof;
	SHPINFO shpIn;
	SHPRECORDHEADER RecordHeader;
    CMapRectangle objRectangle;
    CMapPoint   *pPoint;
	CMapPoints  *pPoints;
	CMapParts   *pParts;
	CMapLine    *pLine;
	CMapPolygon *pPolygon;

	bEof = FALSE;
	switch ( m_shpType )
	{
		case NULLSHP:
			return FALSE;
			break;
        case POINT:
			{

				SHPPTCONTENT point;
				//读入点记录
				for ( i = 1 ; i <= m_iRecordCount ; i++ )
				{
					iLength = SetRecordPos(i);
					if ( iLength <= 0 )
						return FALSE;

					//获得记录头信息
					if (!GetRecordHeader(RecordHeader))
						return FALSE;
					//记录内容长度是否与shp实体大小一致,索引记录长度是否与记录内容长度一致
					if ( RecordHeader.iContentLength*2 != sizeof(point)||
						RecordHeader.iContentLength*2 != iLength)
						return FALSE;
					if(fShp.Read(&point,sizeof(point))!= sizeof(point))
							return FALSE;
					pPoint = new CMapPoint;
					iIndex = i - 1;
				    pPoint->SetIndex(iIndex);
					if ( pPoint == NULL )
						return FALSE;
					//类型是否匹配
					if ( point.iShpType != m_shpType)
						return FALSE;
                    pPoint->SetX(point.dbX );
					pPoint->SetY(point.dbY );
					m_ObList.AddTail(pPoint);

				 }
			}
			break;
        case POLYLINE:

			pszBuffer = new char [MAX_BUFFER_SIZE]; //分配缓冲区
			if ( pszBuffer == NULL )
				return FALSE;
			memset(pszBuffer , 0 , MAX_BUFFER_SIZE);

			//读入线记录
			for ( i = 1 ; i <= m_iRecordCount ; i++ ) //m_iRecordCount
			{

				iLength = SetRecordPos(i);
				if ( iLength <= 0 )
					return FALSE;

				pLine = new CMapLine();
				iIndex = i - 1;
				pLine->SetIndex(iIndex);
				if ( pLine == NULL )
					return FALSE;
				//获得记录头信息
				if (!GetRecordHeader(RecordHeader))
					return FALSE;
				if ( !GetShpInfo(shpIn))
					return FALSE;
				 if (shpIn.ishpType != POLYLINE )//类型不匹配
                 {
                      delete pLine;
					  continue;
				 }
				 if ( shpIn.iNumParts*sizeof(int) > MAX_BUFFER_SIZE )
                 {
					 //多义线段数大于最大缓冲区长度,忽略该对象
					 delete pLine;
					 continue;
                 }
				 //计算对象内容实际长度
				 j = sizeof(SHPINFO) + shpIn.iNumParts*sizeof(int) ;
				 j += shpIn.iNumPoints*sizeof(SHPPOINT);

				 //判断实际长度是否与索引文件中记录的一致
				 if ( RecordHeader.iContentLength*2 != j )
                 {
				     delete pLine;
					 continue;
				 }
				 //设置shp矩形范围
				 objRectangle.SetLeft(shpIn.Box[0].dbX);
				 objRectangle.SetTop(shpIn.Box[0].dbY);
				 objRectangle.SetRight(shpIn.Box[1].dbX);
				 objRectangle.SetBottom(shpIn.Box[1].dbY);
				 pLine->SetExtent(objRectangle);
				 pIParts = new int[shpIn.iNumParts];
				 if ( pIParts == NULL )
                 {
					delete pLine;
					return FALSE;
				 }
				 //读入多义线段索引
				 if ( fShp.Read(pIParts,shpIn.iNumParts*4) != (uint)(shpIn.iNumParts*4))
				 {
					delete pLine;
					return FALSE;
				 }

				 //点坐标存储所占字节数
				 iLength = shpIn.iNumPoints*sizeof(SHPPOINT);
				 //初始化缓冲区数据
				 iDataLen = ReadPoint(pszBuffer,iLength,bEof);
				 if ( iDataLen < 0 )
                 {
					  delete pLine;
					  delete pIParts;
					  return FALSE;

				 }
				 lpVer = pszBuffer;
				 for ( j = 0 ;  j < shpIn.iNumParts ; j++ )
				 {
					pParts  = new CMapParts();
					pPoints = new CMapPoints();
 					if ( pParts == NULL || pPoints == NULL)
						return FALSE;
					if ( j == shpIn.iNumParts - 1 )
                    {
					    k = pIParts[j];        //本段第一个顶点索引
					    m = shpIn.iNumPoints ; //下一个段第一个顶点索引
                    }
              	    else
                    {
						k = pIParts[j];
						m = pIParts[j+1];
					}
					//处理第i段的顶点
				    for ( ; k < m ; k++)
                    {
						pPoint = new CMapPoint();
						if ( pPoint == NULL )
							return FALSE;

						//需要读入数据更新缓冲区
						if ( lpVer == pszBuffer + iDataLen && !bEof)
						{
							iDataLen = ReadPoint(pszBuffer,iLength,bEof);
							if ( iDataLen < 0 )
							{
								delete pPoint;
								delete pPoints;
								delete pLine;
								delete pIParts;

								return FALSE;

							}
							lpVer = pszBuffer;
						}
						dbTmp = *(double*)lpVer;
					    pPoint->SetX(dbTmp);
                        lpVer += 8;
						//需要读入数据更新缓冲区
						if ( lpVer == pszBuffer + iDataLen && !bEof)
						{
							iDataLen = ReadPoint(pszBuffer,iLength,bEof);
							if ( iDataLen < 0 )
							{
								delete pPoint;
								delete pPoints;
								delete pLine;
								delete pIParts;
								return FALSE;

							}
							lpVer = pszBuffer;
						}
						dbTmp = *(double*)(lpVer);
						lpVer += 8;
                        pPoint->SetY(dbTmp);
						pPoints->Add(pPoint);

					}
                    pParts->Add(pPoints);
				    pLine->Add(pParts);

				 }
				 m_ObList.AddTail( pLine);
				 delete []pIParts;

            }
			delete []pszBuffer;

			break;
        case POLYGON:
			pszBuffer = new char [MAX_BUFFER_SIZE]; //分配缓冲区
			if ( pszBuffer == NULL )
				return FALSE;
			memset(pszBuffer , 0 , MAX_BUFFER_SIZE);

			//读入多边形记录
			for ( i = 1 ; i <= m_iRecordCount ; i++ ) //m_iRecordCount
			{
				iLength = SetRecordPos(i);
				if ( iLength <= 0 )
					return FALSE;

				pPolygon = new CMapPolygon();
				iIndex = i - 1;
                pPolygon->SetIndex(iIndex);
				if (pPolygon == NULL )
					return FALSE;
				//获得记录头信息
				if (!GetRecordHeader(RecordHeader))
					return FALSE;
				if ( !GetShpInfo(shpIn))
					return FALSE;
				 if (shpIn.ishpType != POLYGON )//类型不匹配
                 {
                      delete pPolygon;
					  continue;
				 }
				 if ( shpIn.iNumParts*sizeof(int) > MAX_BUFFER_SIZE )
                 {
					 //复合多边型中的多边形个数大于最大缓冲区长度,忽略该对象
					 delete pPolygon;
					 continue;
                 }
				 //计算对象内容实际长度
				 j = sizeof(SHPINFO) + shpIn.iNumParts*sizeof(int) ;
				 j += shpIn.iNumPoints*sizeof(SHPPOINT);

				 //判断实际长度是否与索引文件中记录的一致
				 if ( RecordHeader.iContentLength*2 != j )
                 {
				     delete pPolygon;
					 continue;
				 }
				 //设置shp矩形范围
				 objRectangle.SetLeft(shpIn.Box[0].dbX);
				 objRectangle.SetTop(shpIn.Box[0].dbY);
				 objRectangle.SetRight(shpIn.Box[1].dbX);
				 objRectangle.SetBottom(shpIn.Box[1].dbY);
				 pPolygon->SetExtent(objRectangle);
				 pIParts = new int[shpIn.iNumParts];
				 if ( pIParts == NULL )
                 {
					delete pPolygon;
					return FALSE;
				 }
				 //读入复合多边型段索引
				 if ( fShp.Read(pIParts,shpIn.iNumParts*4) != (uint)(shpIn.iNumParts*4))
				 {
					delete pPolygon;
					return FALSE;
				 }

				 //点坐标存储所占字节数
				 iLength = shpIn.iNumPoints*sizeof(SHPPOINT);
				 //初始化缓冲区数据
				 iDataLen = ReadPoint(pszBuffer,iLength,bEof);
				 if ( iDataLen < 0 )
                 {
					  delete pPolygon;
					  delete pIParts;
					  return FALSE;

				 }
				 lpVer = pszBuffer;
				 for ( j = 0 ;  j < shpIn.iNumParts ; j++ )
				 {
					pParts  = new CMapParts();
					pPoints = new CMapPoints();
 					if ( pParts == NULL || pPoints == NULL)
						return FALSE;
					if ( j == shpIn.iNumParts - 1 )
                    {
					    k = pIParts[j];        //本段第一个顶点索引
					    m = shpIn.iNumPoints ; //下一个段第一个顶点索引
                    }
              	    else
                    {
						k = pIParts[j];
						m = pIParts[j+1];
					}
					//处理第i段的顶点
				    for ( ; k < m ; k++)
                    {
						pPoint = new CMapPoint();
						if ( pPoint == NULL )
							return FALSE;

						//需要读入数据更新缓冲区
						if ( lpVer == pszBuffer + iDataLen && !bEof)
						{
							iDataLen = ReadPoint(pszBuffer,iLength,bEof);
							if ( iDataLen < 0 )
							{
								delete pPolygon;
								delete pIParts;
								return FALSE;

							}
							lpVer = pszBuffer;
						}
						dbTmp = *(double*)lpVer;
					    pPoint->SetX(dbTmp);
                        lpVer += 8;
						//需要读入数据更新缓冲区
						if ( lpVer == pszBuffer + iDataLen && !bEof)
						{
							iDataLen = ReadPoint(pszBuffer,iLength,bEof);
							if ( iDataLen < 0 )
							{
								delete pPolygon;
								delete pIParts;
								return FALSE;

							}
							lpVer = pszBuffer;
						}
						dbTmp = *(double*)(lpVer);
                        pPoint->SetY(dbTmp);
						pPoints->Add(pPoint);
						lpVer += 8;

					}
                    pParts->Add(pPoints);
				    pPolygon->Add(pParts);

				 }
				 m_ObList.AddTail( pPolygon);
			     delete []pIParts;
            }
			delete []pszBuffer;
			break;
		default:
			return FALSE;
			break;

    }
	return TRUE;

}

/*************************************************
  描述:         计算每条shp对象相对文件头的偏移量
  输入:         记录索引值(从零开始)
  输出：        该shp对象数据在文件中的位置
*************************************************/
int CShpFile::SetRecordPos( int iRecord )
{
    unsigned int iOffset,iTmp;
	SHXRECORD shxRD;

	if ( iRecord < 0 )
		return 0;
    //获得索引文件记录偏移量相对文件头
	if (iRecord == 1 )
    	iOffset = sizeof(SHPHEADER)  ;
	else
		iOffset = sizeof(SHPHEADER) + (iRecord-1)*sizeof(shxRecord) ;
	if ( iOffset > m_shxFileLength*2 - sizeof(shxRecord) )
		return 0;
	fShx.Seek( iOffset , CFile::begin );
	int m = sizeof(shxRD);
	fShx.Read( &shxRD , sizeof(shxRD));
	iTmp = shxRD.iOffset;
	SwapWord(sizeof(int),&iTmp);
	fShp.Seek(iTmp*2 ,  CFile::begin );
	iTmp = shxRD.iContentLength;
	SwapWord(sizeof(int),&iTmp);


	return iTmp*2;
}

/*************************************************
  描述:         获得每条shp对象记录记录头的信息
  输入:         记录头结构对象
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/

BOOL CShpFile::GetRecordHeader(SHPRECORDHEADER& RecordHeader )
{
	int iLength,iNum;

	if(fShp.Read(&RecordHeader,sizeof(RecordHeader))!= sizeof(RecordHeader))
		return FALSE;
	if ( !m_bBigEndian )
    {
		iNum    = RecordHeader.iRecordNum;
        iLength = RecordHeader.iContentLength;
		SwapWord(sizeof(int),&iLength);
		SwapWord(sizeof(int),&iNum);
		RecordHeader.iRecordNum = iNum;
        RecordHeader.iContentLength = iLength;
	}
    return TRUE;


}
/*************************************************
  描述:         获得每条shp对象描述信息
  输入:         描述信息结构对象
  输出：        成功返回TRUE 失败返回FALSE
*************************************************/

BOOL CShpFile::GetShpInfo(SHPINFO& varInfo)
{

	if(fShp.Read(&varInfo,sizeof(varInfo))!= sizeof(varInfo))
		return FALSE;
	return TRUE;
}

/*************************************************
  描述:         读入点对象数据
  输入:         数据缓冲区指针 缓冲区最大32K
                如果超出需要分多次读取，要读取的长度、
				是否已读取完成
  输出：        读取数据的实际长度
*************************************************/

int  CShpFile::ReadPoint(char* pszBuffer,int& iLength,BOOL& bEof)
{
	if ( iLength > MAX_BUFFER_SIZE)
    {
		iLength -= MAX_BUFFER_SIZE;
		if ( fShp.Read(pszBuffer,MAX_BUFFER_SIZE) != MAX_BUFFER_SIZE )
			return FILE_READERR;
		bEof = FALSE;
		return MAX_BUFFER_SIZE;
    }
    else
    {
		if ( fShp.Read(pszBuffer,iLength) != (uint)iLength )
			return FILE_READERR;
		bEof = TRUE;
		return iLength;
	}
}


/*************************************************
  描述:         根据点对象查找shp对象是否选中
  输入:         点对象
  输出：        查找到shp对象的索引值 返回值-表示未查找到
*************************************************/
int CShpFile::SearchShape(CMapPoint& pt )
{
	unsigned long iCount;
	POSITION pos;
    CMapPolygon *pPolygon;

	if ( GetShpType() != POLYGON ) //只判断多边形对象
		return -1;
	iCount = m_ObList.GetCount()-1;

	for ( pos = m_ObList.GetHeadPosition()  ; pos != NULL ;  )
    {
		pPolygon = (CMapPolygon*)m_ObList.GetAt(pos);
		if ( pPolygon->IsPointIn(pt) )
		 	return pPolygon->GetIndex();
		m_ObList.GetNext(pos);
    }
	return -1;

}

/*************************************************
  描述:         绘制shp对象
  输入:         
  输出：        osg::Group 指针
*************************************************/
osg::ref_ptr<osg::Group>  CShpFile::DrawShp(int mode)
{           
	osg::ref_ptr<osg::Group> gp = NULL;
	switch ( m_shpType )
	{
		case POINT:
		{
			DrawPoint(mode);	currShpType = 0;
		}
		break;
		
		case POLYGON: 
		{
			gp = DrawPolygon(mode);		currShpType = 2;
		}  
		break;   	 

		case POLYLINE: 
		{
			gp = DrawPLine(mode);		currShpType = 1;
		}  
		break;   	 

		default:
		break;
	}
	return gp;
}

/*************************************************
  描述:         绘制点对象
  输入:         设备指针、图例对象指针、坐标变换参数结构对象
  输出：        无
*************************************************/
std::string DeleteSpace(CString csValue)
{
	/* 如果csValue行尾有空格，则把它清除掉 */
	char namestr[128];

	std::string name = csValue;
	sprintf(namestr, "%s", name.c_str());
	for(unsigned int i=0; i<strlen(namestr); i++)
	{
		if((namestr[i]==' ') || (namestr[i]==0x9))  namestr[i] = 0;
	}
	name = namestr;
	
	return name;
}

void InsertCultPoint(SimPoint &sp)
{
	char tile[16];
	Bucket bt;
	set_bucket(&bt, sp.lon, sp.lat);
	gen_index_str(&bt, tile);
	sp.tile = tile;
	
	/* 检查sp.name，是不是灯点没加Lod的名称 */
	GetWidthFromName(sp.name);
	if((nameParts[0] == "USR") && (nameParts[2] == "LIGHTS"))
	{
		if(nameParts.size() <= 4) sp.name = sp.name + "_0_30000";
		sp.name = sp.name + "_" + sp.tile;
		CultPoints.push_back(sp);
	}
	
	/* 如果是点树，加到TreePoints里面 */
	if((nameParts[0] == "Forest") || (nameParts[0] == "forest") || (nameParts[0] == "EvergreenForest"))
	{
		sp.name = getTreeName(sp.name);
#if 1
		/* 调整程序，把树降低一点高度。降高1.5米 */
		Geodesy geod;
		Carton  cart;
		geod._lon = sp.lon; geod._lat = sp.lat; geod._elevation = sp.elev;
		geod._elevation -= 1.5;
		GeodToCart(&geod, &cart);
		sp.x  = cart.x; sp.y  = cart.y; sp.z  = cart.z; sp.elev = geod._elevation;
#endif
		TreePoints.push_back(sp);
	}	
}


void CShpFile::DrawPoint(int mode)
{
	int    iIndex,iField;
	CMapFields *pFields;
	CMapField  *pField;
	CMapPoint *pPoint;
	POSITION pos;
	CString csValue;
	CString csName;
	Geodesy geod;
	Carton  cart;
	
	for( pos = m_ObList.GetHeadPosition(); pos != NULL; )
	{  
		SimPoint sp;
		pPoint = (CMapPoint*)m_ObList.GetAt(pos);

		iIndex = pPoint->GetIndex();
		pRecordSet->Move(iIndex,BookmarkFirst);
		pFields = pRecordSet->GetFields(0);
		
		/* 取得名字 */
		iField  = 0;
		pField = pFields->GetField(iField);
		csName = pField->GetName();
		while(csName != "NAME")
		{
			iField++;
			pField = pFields->GetField(iField);
			csName = pField->GetName();
		}
		csValue = pField->GetValueAsString();
		sp.name = DeleteSpace(csValue);

		/* 取得X */	/* 取得Y */
		DrawPointElement(pPoint, &sp);

		/* 查找极值 */
		SetMMLL(sp.lon, sp.lat);

		/* 取得高程 */
		if(mode == NEED_ELEVATION)
		{
			geod._lon = sp.lon;
			geod._lat = sp.lat;
			sp.elev = geod._elevation = GetElevation(geod._lon, geod._lat, 0.0) + 5.0;
			GeodToCart(&geod, &cart);
			sp.x = cart.x; sp.y = cart.y; sp.z = cart.z;
			InsertCultPoint(sp);
		}

		m_ObList.GetNext(pos);
	}
}

/*************************************************
  描述:         绘制单个的点对象
  输入:         设备指针、点对象对象指针、
                坐标变换参数结构对象
  输出：        无
*************************************************/
void CShpFile::DrawPointElement(CMapPoint *pPoint, SimPoint *sp)
{
	double iX,iY;

	if ( pPoint == NULL )
		return;
	
	iX = pPoint->GetX();
	iY = pPoint->GetY();
	
	sp->lon = iX;
	sp->lat = iY;
	
}

/*************************************************
  插入细节纹理名称
*************************************************/
void InsertDetailedTxt(std::string dt)
{
	unsigned int i;
	int iFlag = 0;
	
	if(dt.substr(0, 11) == "DETAILEDTXT")
	{
		for(i=0; i<allDetailedTxt.size(); i++)
		{
			if(dt == allDetailedTxt[i]) {	iFlag = 1; break; }
		}
		if(!iFlag)	allDetailedTxt.push_back(dt);
		ifHasIslandInCurrentTerrain = 1;
	}
}

/*************************************************
  描述:         根据多边形对象创建Group
  输入:         多边形对象
  输出：        osg::Group
*************************************************/
osg::ref_ptr<osg::Group> CShpFile::DrawPolygon(int mode)
{
	int    iIndex,iField; 
	CMapFields *pFields;
	CMapField  *pField;
	POSITION pos;
	CString csValue;
	CString csName;
	CMapPolygon *pPolygon;

	osg::ref_ptr<osg::Group> gp = new osg::Group;
	osg::ref_ptr<osg::Group> gd;

	for( pos = m_ObList.GetHeadPosition(); pos != NULL; )
	{  
		pPolygon = (CMapPolygon*)m_ObList.GetAt(pos);
		iIndex 	= pPolygon ->GetIndex();
		pRecordSet->Move(iIndex,BookmarkFirst);
		pFields = pRecordSet->GetFields(0);
		
		/* 首先找到FieldName为 “NAME” 的项*/
		iField  = 0;
		pField = pFields->GetField(iField);
		csName = pField->GetName();
		while(csName != "NAME")
		{
			iField++;
			pField = pFields->GetField(iField);
			csName = pField->GetName();
		}				
		csValue = pField->GetValueAsString();

		/* 判断并插入细节纹理名称 */
		std::string dtname = DeleteSpace(csValue);
		InsertDetailedTxt(dtname);

		/* 如果当前名称是Airport，那么获取机场的高度 */
		isAirport = false;
		if((dtname.substr(0, 7) == DEFINE_AIRPORT_STRING) || (dtname.substr(0, 7) == "AIRPORT"))
		{
			GetWidthFromName(dtname);
			if(nameParts.size() >= 2)
				m_elevOfAirport = atof(nameParts[1].c_str());
			else
				m_elevOfAirport = 0.0;
			isAirport = true;
		}

		/* 如果当前名称是Lake，那么获取湖面的高度 */
		if((dtname.substr(0, 4) == DEFINE_LAKE_STRING) || (dtname.substr(0, 4) == "LAKE"))
		{
			GetWidthFromName(dtname);
			if(nameParts.size() >= 2)
			{
				m_elevOfLake = atof(nameParts[1].c_str()) - 1.9;
				dtname = nameParts[0];
			}
			else
				m_elevOfLake = 0.0;
		}

		
		/* 如果当前Polygon是岛，设置全局变量 */
		if ((csValue[0] == 'I') && (csValue[1] == 's') && (csValue[2] == 'l') && (csValue[3] == 'a') && (csValue[4] == 'n') && (csValue[5] == 'd')) 
			ifHasIslandInCurrentTerrain = 1;
		
		/* 看当前Polygon是否是湖 */
		if(	((csValue[0] == 'L') && (csValue[1] == 'a') && (csValue[2] == 'k') && (csValue[3] == 'e')) ||
				((csValue[0] == 'I') && (csValue[1] == 's') && (csValue[2] == 'l') && (csValue[3] == 'a') && (csValue[4] == 'n') && (csValue[5] == 'd')) )
			isLake = true;
		if((csValue[0] == 'I') && (csValue[1] == 's') && (csValue[2] == 'l') && (csValue[3] == 'a') && (csValue[4] == 'n') && (csValue[5] == 'd') && (csValue[6] == 'I') && (csValue[7] == 'n') && (csValue[8] == 'O') && (csValue[9] == 'c') && (csValue[10] == 'e') && (csValue[11] == 'a') ) 
		{
			isLake = false;
			islandInOcean = true;
		}
		
		gd = DrawPolygonElement(mode, pPolygon);
		for(unsigned int i = 0; i<gd->getNumChildren(); i++)
		{
			osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gd->getChild(i));
			gnode->setName(dtname);
			gp->addChild(gnode.get());
		}

		m_ObList.GetNext(pos);
	}
	
	return gp.release();
}


/*************************************************
  描述:         根据多边形对象创建Geode
  输入:         多边形对象
  输出：        osg::Geode
*************************************************/
void CShpFile::GetCurrPolygonPointElevation(Geodesy *geod, int mode)
{
	if(mode == NEED_ELEVATION)
	{
		if(islandInOcean)
			geod->_elevation = 0.0;
		else if((isLake) && (m_elevOfLake == 0.0))
			m_elevOfLake = geod->_elevation = GetElevation(geod->_lon, geod->_lat, 0.0);
		else if((isLake) && (m_elevOfLake != 0.0))
			geod->_elevation = m_elevOfLake;
		else if((isAirport) && (m_elevOfAirport != 0.0))
			geod->_elevation = m_elevOfAirport;
		else if(!isLake)
			geod->_elevation = GetElevation(geod->_lon, geod->_lat, 0.0);
	}
	else
		geod->_elevation = 0.0;
}

osg::ref_ptr<osg::Group> CShpFile::DrawPolygonElement(int mode, CMapPolygon* pPolygon)
{
	int i,j,k, isLastPointValid;
	CMapPoint *pPoint;
	CMapPoints *pPoints;
	CMapParts *pParts;
	Geodesy geod;
	Carton  cart;
	osg::Vec3 lastPoint;
	osg::Vec3 currPoint; 

	if ( pPolygon == NULL )	return NULL;

	osg::ref_ptr<osg::Group> gp = new osg::Group;
	for (i = 0 ; i < pPolygon->GetCount() ; i++ )
	{
		osg::ref_ptr<osg::Geode> gd = new osg::Geode;
		gp->addChild(gd.get());
		gd->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
		osg::ref_ptr<osg::Geometry>	gm = new osg::Geometry;
		osg::ref_ptr<osg::Vec3Array> va = new osg::Vec3Array;
		gm->setVertexArray(va);
		isLastPointValid = 0;

		pParts =  pPolygon->GetParts(i);
		for ( j = 0 ; j < pParts->GetCount() ; j++)
		{
			pPoints = pParts->Get(j);
			for ( k = 0 ; k < pPoints->GetCount(); k++ )
			{
				pPoint = pPoints->Get(k);
				geod._lon = pPoint->GetX();
				geod._lat = pPoint->GetY();
				SetMMLL(geod._lon, geod._lat);
				
				currPoint.set(geod._lon, geod._lat, 0.0);
				if(isLastPointValid)
				{
					/* 首先在上一点和当前点之间按照固定的距离插入若干中间点。*/
					std::vector<osg::Vec3> vec_ab;
					interpolation(lastPoint, currPoint, vec_ab, INSERTDISTANCE);
					vec_ab.push_back(currPoint);
					
					/* 把每个插入点都计算高程，并插入点队列中 */
					for(unsigned int m=0; m<vec_ab.size(); m++)
					{
						geod._lon = vec_ab[m].x(); geod._lat = vec_ab[m].y();
						GetCurrPolygonPointElevation(&geod, mode);
						GeodToCart(&geod, &cart);
	
						osg::Vec3 v; 
						v.set(cart.x, cart.y, cart.z);
						va->push_back(v);
					}
				}else
				{
					GetCurrPolygonPointElevation(&geod, mode);
					GeodToCart(&geod, &cart);

					osg::Vec3 v; 
					v.set(cart.x, cart.y, cart.z);
					va->push_back(v);
				}

				/* 插入点的起始点 */
				lastPoint = currPoint;
				isLastPointValid = USE_INSERT_VERTEX;
			}
		}
	
		if(va->size() != 0)
		{
			gm->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, va->size()) );
			gd->addDrawable(gm.get());
		}
	}

#if 0
	sprintf(debugFile, "%s\\PolygonElement%d.ive", tempDirectory, elemNum++);
	osgDB::writeNodeFile(*gp, debugFile);
#endif

	return gp.release();
}

/*************************************************
  描述:         解析线的名称，以下划线或者中线作为分段标志，将名称分为几段，
  输入:         string name
  输出：        int width
*************************************************/
float GetWidthFromName(std::string matname)
{
	char matpart[128], mat[256];
	int j = 0, iWidth;
	
	/* 把matname以下划线作为分隔符，分成若干小段，保存在nameParts里面 */
	sprintf(mat, "%s", matname.c_str());
	nameParts.clear();
	for(unsigned int i = 0; i<strlen(mat); i++)
	{
		if((mat[i] != '_')&&(mat[i] != '-')) { matpart[j++] = mat[i];	}
		if((mat[i] == '_')||(mat[i] == '-')) { matpart[j] = 0;	nameParts.push_back(matpart);	j = 0; }
	}
	matpart[j] = 0;	nameParts.push_back(matpart);	j = 0;
	
	/* 根据nameParts的第二项作为宽度值 */ 
	if(nameParts.size() >= 2)
	{
		iWidth = atof(nameParts[1].c_str());
		if((iWidth <= 0) || (iWidth > 100))	iWidth = mRoadWidthInShp;
	}else{
		iWidth = mRoadWidthInShp;
	}
	
	return iWidth;
}


/*************************************************
  描述:         根据多边形对象创建Group
  输入:         多边形对象
  输出：        osg::Group
*************************************************/
osg::ref_ptr<osg::Group> CShpFile::DrawPLine(int mode)
{
	int    iIndex,iField; 
	CMapFields *pFields;
	CMapField  *pField;
	POSITION pos;
	CString csValue;
	CString csName;
	CMapLine *pPLine;
	int  curType = 0;			/* 0:Road; 1:Tree */

	osg::ref_ptr<osg::Group> gp = new osg::Group;
	osg::ref_ptr<osg::Group> gd;

	for( pos = m_ObList.GetHeadPosition(); pos != NULL; )
	{  
		pPLine = (CMapLine*)m_ObList.GetAt(pos);
		iIndex 	= pPLine ->GetIndex();
		pRecordSet->Move(iIndex,BookmarkFirst);
		pFields = pRecordSet->GetFields(0);

		/* 首先找到FieldName为 “NAME” 的项*/
		iField  = 0;
		pField = pFields->GetField(iField);
		csName = pField->GetName();
		while(csName != "NAME")
		{
			iField++;
			pField = pFields->GetField(iField);
			csName = pField->GetName();
		}
		csValue = pField->GetValueAsString();
		
		/* 解析这个名字，从中得到要做的路的宽度值，介于0~100之间，规定名字为Road_???，???为宽度值  */
		std::string dtname = DeleteSpace(csValue);
		float curWidth = GetWidthFromName(dtname);
		
		/* 根据名字来判断，当前做的是路还是线树 */
		if((dtname.substr(0, 4) == "Road")||(dtname.substr(0, 4) == "road"))
			gd = DrawPLineElement(mode, pPLine, curWidth, 0);
		if((dtname.substr(0, 6) == "Forest") || (dtname.substr(0, 6) == "forest"))
			gd = DrawPLineElement(mode, pPLine, curWidth, 1);
		
		if(gd != NULL)
		{
			for(unsigned int i = 0; i<gd->getNumChildren(); i++)
			{
				osg::ref_ptr<osg::Geode> gnode = dynamic_cast<osg::Geode*> (gd->getChild(i));
				gp->addChild(gnode.get());
			}
		}

		m_ObList.GetNext(pos);
	}
	
	return gp.release();
}


/*************************************************
  描述:         根据PolyLine对象创建Geode
  输入:         多边形对象
  输出：        osg::Geode
*************************************************/

/* 已知一条线段上的A, B两点，和需要将AB线段分段的段数，求：靠近A点的最近的一个分段点。*/
osg::Vec3d GetOffsetVertex(osg::Vec3d A, osg::Vec3d B, float distance)
{
	double delx, dely, delz, step, delta;
	int count;
	osg::Vec3d iVec;

	delx = B[0]-A[0];
	dely = B[1]-A[1];
	delz = B[2]-A[2];
	step = distance;

	if (fabs(delx)>fabs(dely))
	{
		count = fabs(delx) / step + 1;
		delta = delx / count;

		iVec[0] = A[0] + delta;
		iVec[1] = (delta)*dely/delx + A[1];
		iVec[2] = (delta)*delz/delx + A[2];
	}
	else if (fabs(dely)>fabs(delx))
	{
		count = fabs(dely) / step + 1;
		delta = dely / count;

		iVec[1] = A[1] + delta;
		iVec[0] = (delta)*delx/dely + A[0];
		iVec[2] = (delta)*delz/dely + A[2];
	}
	else if((fabs(delx) == 0.0) && (fabs(dely) == 0.0))
	{
		iVec[0] = A[0];
		iVec[1] = A[1];
		iVec[2] = A[2];
	}
	else
	{
		/* fabs(delx)和fabs(dely)相等，但都不等于0.0 */
		count = fabs(delx) / step + 1;
		delta = delx / count;

		iVec[0] = A[0] + delta;
		iVec[1] = (delta)*dely/delx + A[1];
		iVec[2] = (delta)*delz/delx + A[2];
	}
	
	return iVec;
}

/* 求取空间两个点A和B之间在直线AB上的的若干(segs)个间隔相等的点。 */
void GetLightPots(osg::Vec3d A, osg::Vec3d B, float distance)
{
	double delx, dely, delz, step, delta;
	std::vector<osg::Vec3d> vec_list;
	int count;
	Geodesy geod;
	Carton  cart;
	
	/* 加入A点 */
	SimPoint sp;
	sp.name = "USR_WHITE_LIGHTS_0_30000";
	cart.x = sp.x = A.x();	cart.y = sp.y = A.y();	cart.z = sp.z = A.z();
	CartToGeod(&cart, &geod);
	sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
	InsertCultPoint(sp);
	
	/* AB之间均分segs份，间隔点都插入到vec_list里面 */
	delx = B[0]-A[0];
	dely = B[1]-A[1];
	delz = B[2]-A[2];
	step = distance;

	if (fabs(delx)>fabs(dely))
	{
		count = fabs(delx) / step + 1;
		delta = delx / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3d iVec;
			iVec[0] = A[0] + delta*(i+1);
			iVec[1] = (delta*(i+1))*dely/delx + A[1];
			iVec[2] = (delta*(i+1))*delz/delx + A[2];
			vec_list.push_back(iVec);
		}
	}
	else if (fabs(dely)>fabs(delx))
	{
		count = fabs(dely) / step + 1;
		delta = dely / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3d iVec;
			iVec[1] = A[1] + delta*(i+1);
			iVec[0] = (delta*(i+1))*delx/dely + A[0];
			iVec[2] = (delta*(i+1))*delz/dely + A[2];
			vec_list.push_back(iVec);
		}
	}
	else if((fabs(delx) == 0.0) && (fabs(dely) == 0.0))
	{
		osg::Vec3d iVec;
		iVec[0] = A[0];
		iVec[1] = A[1];
		iVec[2] = A[2];
		vec_list.push_back(iVec);
	}
	else
	{
		/* fabs(delx)和fabs(dely)相等，但都不等于0.0 */
		count = fabs(delx) / step + 1;
		delta = delx / count;
		for (int i=0;i<count-1;i++)
		{
			osg::Vec3d iVec;
			iVec[0] = A[0] + delta*(i+1);
			iVec[1] = (delta*(i+1))*dely/delx + A[1];
			iVec[2] = (delta*(i+1))*delz/delx + A[2];
			vec_list.push_back(iVec);
		}
	}
	
	/* 把vec_list里面的所有点都加入CultPoint */
	for(unsigned int i=0; i<vec_list.size(); i++)
	{
		osg::Vec3d C = vec_list[i];
		
		sp.name = "USR_WHITE_LIGHTS_0_30000";
		cart.x = sp.x = C.x();	cart.y = sp.y = C.y();	cart.z = sp.z = C.z();
		CartToGeod(&cart, &geod);
		sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
		InsertCultPoint(sp);
	}
	
	/* 最后加入B点 */
	sp.name = "USR_WHITE_LIGHTS_0_30000";
	cart.x = sp.x = B.x();	cart.y = sp.y = B.y();	cart.z = sp.z = B.z();
	CartToGeod(&cart, &geod);
	sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
	InsertCultPoint(sp);	
}

/* 在当前边上插入灯点 */
double totalOffset = 0.0;
void InsertLightPotsOnCurrEdge(std::vector<GeoEdge>::iterator v, float distance)
{
	double currDist  = v->dist;							//当前这条边的长度。
	double validDist;
	int	   numLights;
	osg::Vec3d A = v->Ac, B = v->Bc;

	if(totalOffset > distance / 10.0)
	{
		/* 一、计算出v线段上靠近A端，与B端距离为(distance - totalOffset)的点的值: midA */
		osg::Vec3d midA = GetOffsetVertex(A, B, (distance - totalOffset));
			
		/* 二、计算出midA与B之间能插入多少灯点 */
		validDist = currDist - (distance - totalOffset);
		numLights = (int)( validDist / distance);
		
		/* 三、计算插入numLights点之后，剩下多少距离 */
		double currOffset = validDist - numLights * distance;
		if(currOffset > distance / 10.0)
		{
			/* 先计算出v线段上靠近B端，与A端距离为currOffset的点的值: midC */
			osg::Vec3d midB = GetOffsetVertex(B, A, currOffset);
				
			/* 再在线段A<->midC之间均分numLights份，就是路灯点了，并转存到CultPoints*/
			GetLightPots(midA, midB, distance);
			
			/* 剩下的一段转存为totalOffset，供下一条线段使用 */
			totalOffset = currOffset;
		}else{
			/* currOffset 忽略不计 */
			/* 直接在线段A<->B之间均分numLights份，就是路灯点了，并转存到CultPoints*/
			GetLightPots(midA, B, distance);
		}
	}else{
		/* 忽略不计 totalOffset */
		/* 一、计算出A与B之间能插入多少灯点 */
		validDist = currDist;
		numLights = (int)( validDist / distance);
		
		/* 二、计算插入numLights点之后，剩下多少距离 */
		double currOffset = validDist - numLights * distance;
		if(currOffset > distance / 10.0)
		{
			/* 先计算出v线段上靠近B端，与A端距离为currOffset的点的值: midC */
			int  numOffsetSeg = (int)(validDist / currOffset);
			osg::Vec3d midB = GetOffsetVertex(B, A, numOffsetSeg);
				
			/* 再在线段A<->midC之间均分numLights份，就是路灯点了，并转存到CultPoints*/
			GetLightPots(A, midB, distance);
			
			/* 剩下的一段转存为totalOffset，供下一条线段使用 */
			totalOffset = currOffset;
		}else{
			/* currOffset 忽略不计 */
			/* 直接在线段A<->B之间均分numLights份，就是路灯点了，并转存到CultPoints*/
			GetLightPots(A, B, distance);
		}
	}
}

/* 在边线上以一定的间隔生成灯点 */
int CreateLightPotOnSideLine(float distance)
{
	totalOffset = 0.0;
	/* 根据经纬高格式的边计算插入点 */
	std::vector<GeoEdge>::iterator v = SideGeoLine.begin();
	while(v != SideGeoLine.end())
	{
		double currDist  = v->dist;							//当前这条边的长度。
		int	   numLights = (int)( (currDist + totalOffset) / distance);		//当前这条边应当有的灯点数
		if(numLights > 0)
		{
			InsertLightPotsOnCurrEdge(v, distance);
		}else{
			totalOffset += currDist;
		}
		
		/* 下一条边 */
		v++;
	}
	
	return 1;	
}

/* 生成一条绿化带 */
osg::ref_ptr<osg::Group> CreateTreeLine(float width)
{
	/* 对于所有的线做循环 */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	for(std::vector<MyPolyLine>::iterator tt = totallines.begin(); tt != totallines.end(); tt++)
	{
		/* 左边 */
		/* 在已知线左边做成一条绿化带 */
		MyPolyLine mpl = *tt;
		osg::ref_ptr<osg::Geode> gLTree = CreateRoadWithPolyLine(&mpl, width/2.0, 0, false);
		if(gLTree.valid())
		{
			gLTree->setName(DEFINE_FOREST_STRING);
			gp->addChild(gLTree.get());
		}

		/* 右边 */
		/* 在已知线右边做一条绿化带 */
		osg::ref_ptr<osg::Geode> gRTree = CreateRoadWithPolyLine(&mpl, width/2.0, 1, false);
		if(gRTree.valid())
		{
			gRTree->setName(DEFINE_FOREST_STRING);
			gp->addChild(gRTree.get());
		}
	}
	
	return gp.release();
}

/* 生成一条两边带绿化带，并且有路灯的路的文化信息 */
osg::ref_ptr<osg::Group> CreateRoadsWithTree(float width)
{
	/* 对于所有的线做循环 */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	for(std::vector<MyPolyLine>::iterator tt = totallines.begin(); tt != totallines.end(); tt++)
	{
		/* 左边 */
		/* 先把已知线做成一条路 */
		MyPolyLine mpl = *tt;
		osg::ref_ptr<osg::Geode> gLRoad = CreateRoadWithPolyLine(&mpl, width/2.0, 0, true);
		if(gLRoad.valid())
		{
#if 0
			osgDB::writeNodeFile(*gLRoad, "F:\\glroad.ive");
			gLRoad->setName(DEFINE_ROAD_STRING);
			gp->addChild(gLRoad.get());
#endif
		}
		
		/* 根据这个Geode找到平行那根线，即路对面的边，保存在SideLine里面 */
		if(gLRoad.valid())	CreateSideLineWithCurrentRoad(gLRoad);
		
#if 0
	/* 暂不生成道路上的随机灯点。。2013-01-22 */
		/* 根据路对面的边线SideLine，创建相等间隔的灯点 */
		if((mLightPotDistance > 0.0) && (CurrentLampType == true))
		{
			if(CreateLightPotOnSideLine(mLightPotDistance))
			{
			}
		}
	/* 注释完毕 */	
#endif	//0
					
#if 0
	/* 暂不生成道路边的绿化带。。2013-01-23 */
		/* 根据对边线做出左边绿化带 */
		osg::ref_ptr<osg::Geode> gLeft  = CreateRoadWithPolyLine(&SideLine, 20.0, 0, false);
		if(gLeft.valid())
		{
			gLeft->setName(DEFINE_FOREST_STRING);
			gp->addChild(gLeft.get());
		}
	/* 注释完毕 */	
#endif	//0

		/* 右边 */
		/* 先把已知线做成一条路 */
		osg::ref_ptr<osg::Geode> gRRoad = CreateRoadWithPolyLine(&mpl, width/2.0, 1, true);
		if(gRRoad.valid())
		{
			//osgDB::writeNodeFile(*gLRoad, "F:\\grroad.osg");
			//gRRoad->setName(DEFINE_ROAD_STRING);
			//gp->addChild(gRRoad.get());
		}
		
		/* 根据这个Geode找到平行那根线，即路对面的边，保存在SideLine里面 */
		if(gRRoad.valid())	CreateSideLineWithCurrentRoad(gRRoad);

#if 0
	/* 暂不生成道路上的随机灯点。。2013-01-22 */
		/* 根据路对面的边线SideLine，创建相等间隔的灯点 */
		if((mLightPotDistance > 0.0) && (CurrentLampType == true))
		{
			if(CreateLightPotOnSideLine(mLightPotDistance))
			{
			}
		}
	/* 注释完毕 */	
#endif	//0
				
#if 0
	/* 暂不生成道路边的绿化带。。2013-01-23 */
		/* 根据对边线做出右边绿化带 */
		osg::ref_ptr<osg::Geode> gRight  = CreateRoadWithPolyLine(&SideLine, 20.0, 1, false);
		if(gRight.valid())
		{
			gRight->setName(DEFINE_FOREST_STRING);
			gp->addChild(gRight.get());
		}
	/* 注释完毕 */	
#endif	//0
		
		/* 合并左边路和右边路，组合成一条大路 */
		if((gLRoad.valid()) && (gRRoad.valid()))
		{
			osg::ref_ptr<osg::Geode> gRoad = MergeRoadFromLAndR(gLRoad, gRRoad, 1);
#if 0
			osgDB::writeNodeFile(*gRoad, "F:\\groad.ive");
#endif
			gRoad->setName(DEFINE_ROAD_STRING);
			gp->addChild(gRoad.get());
		}
	}
	
	return gp.release();
}

/* 根据当前的路的中线，生成路的边线 */
void CreateSideLineWithCurrentRoad(osg::ref_ptr<osg::Geode> gd)
{
	Geodesy geodA, geodB, gSL;
	Carton  cartA, cartB, cSL;
	osg::Vec3 st, ed;

	SideLine.clear();	SideGeoLine.clear();
	osg::Geometry *g = dynamic_cast<osg::Geometry *> (gd->getDrawable(0));
	osg::Vec3Array *va = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	
	/* 把里面奇数点提出来，弄成一个MyPolyLine */
	osg::Vec3Array::iterator v = va->begin();	v++;
	osg::Vec3 start; start.set(*v);	 			v+=2;
	osg::Vec3 startGeo; startGeo.set(start);

	/*当前点转为海拔高度 */
	cSL.x = start.x();	cSL.y = start.y(); cSL.z = start.z();	CartToGeod(&cSL, &gSL);	gSL._elevation = 0.0;
	GeodToCart(&gSL, &cSL);	start.x() = cSL.x;	start.y() = cSL.y;	start.z() = cSL.z;

	int i=1;
	for(v; v != va->end(); v++)
	{
		if(i%2)
		{
			/* 计算直角坐标系的边 */
			osg::Vec3 end; 	end.set(*v);
			osg::Vec3 endGeo; endGeo.set(end);

			/*当前点转为海拔高度 */
			cSL.x = end.x();	cSL.y = end.y(); cSL.z = end.z();	CartToGeod(&cSL, &gSL);	gSL._elevation = 0.0;
			GeodToCart(&gSL, &cSL);	end.x() = cSL.x;	end.y() = cSL.y;	end.z() = cSL.z;

			SimEdge se;
			se.A = start; se.B = end; se.x = 2*i - 1; se.y = 2*i + 1; se.z = 0;
			SideLine.push_back(se);	
			
			/* 计算经纬高坐标系的边 */
			cartA.x = startGeo.x();  	cartA.y = startGeo.y();  	cartA.z = startGeo.z();
			CartToGeod(&cartA, &geodA);
			cartB.x = endGeo.x();  	cartB.y = endGeo.y();  	cartB.z = endGeo.z();
			CartToGeod(&cartB, &geodB);
			GeoEdge ge; ge.A = geodA; ge.B = geodB; ge.Ac = start; ge.Bc = end;	ge.distance();
			SideGeoLine.push_back(ge);

			/* 下一条边 */
			start = end;	startGeo = endGeo;
		}
		i++;
	}

	return;
}

/* 根据当前的路的中线，生成一条路的多边形序列，即路面，“根据路中线生成路面” */
/* 这里的 line，应当是海拔高度为0.0的线。而返回的Geode, 则是计算了海拔高度的线。 */
osg::ref_ptr<osg::Geode> CreateRoadWithPolyLine(MyPolyLine *line, float width, int Direction, bool hasLight)
{
	Geodesy geod;
	Carton  cart;

	if(line->size() == 0) return NULL;

	/* 根据已知的一条线，计算出一组相连的矩形，这组矩形相同一侧的边是这条已知线，且这组矩形的宽度相等。*/
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	for(MyPolyLine::iterator p = line->begin(); p != line->end(); p++)
	{
		SimEdge se = *p;
		osg::ref_ptr<osg::Geometry> gm;
		gm = CreateRectangleWithEdge(&se, width, Direction);
		gd->addDrawable(gm.get());		

	}

#if 0
	osgDB::writeNodeFile(*gd, "F:\\g1.ive");
#endif


	/* 求相邻两个矩形之间的公共边 */
	osg::Geometry *g       = dynamic_cast<osg::Geometry *> (gd->getDrawable(0));
	osg::Vec3Array *coords = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	osg::Vec3 A; A.set((*coords)[2]);
	osg::Vec3 B; B.set((*coords)[3]);

	for(unsigned int ii = 1; ii < gd->getNumDrawables(); ii++)
	{
		osg::Geometry *gy       = dynamic_cast<osg::Geometry *> (gd->getDrawable(ii));
		osg::Vec3Array *coordd  = dynamic_cast<osg::Vec3Array *> (gy->getVertexArray());
		osg::Vec3 C; C.set((*coordd)[2]);
		osg::Vec3 D; D.set((*coordd)[3]);
		
		osg::Vec3 I;
		double dist = LenOf2V(B, C);
		if(dist < 1.0) 
			I = B;
		else 
			intersect(A, B, C, D, I);
		(*coords)[3].set(I); (*coordd)[2].set(I);
		
		A = C; B = D; coords = coordd;
	}

	/* 组合成一个紧凑的Geode，并增加纹理坐标，纹理坐标增加方法是每个四边形用一个整图。 */
	osg::ref_ptr<osg::Geometry> gm2 = new osg::Geometry;
	osg::ref_ptr<osg::Vec3Array> va = new osg::Vec3Array;
	gm2->setVertexArray(va);
#if 0
	osg::ref_ptr<osg::Vec2Array> ta = new osg::Vec2Array;
	gm2->setTexCoordArray(0, ta);
	osg::Vec2 tStart[2], tEnd[2];	tStart[0].set(0.0, 0.0);	tStart[1].set(1.0, 0.0);	tEnd[0].set(0.0, 1.0);		tEnd[1].set(1.0, 1.0);
#endif
	osg::Vec3Array *coordd;
	int iFlag = 0;
	for(unsigned int ii = 0; ii < gd->getNumDrawables(); ii++)
	{
		osg::Geometry *gy       = dynamic_cast<osg::Geometry *> (gd->getDrawable(ii));
		coordd  = dynamic_cast<osg::Vec3Array *> (gy->getVertexArray());

		va->push_back((*coordd)[0]);
		va->push_back((*coordd)[2]);
		
#if 0
		if(iFlag == 0)	{	ta->push_back(tStart[0]);	ta->push_back(tStart[1]);	}
		if(iFlag == 1)	{	ta->push_back(tEnd[0]);		ta->push_back(tEnd[1]);		}
		if(iFlag == 0)	iFlag = 1; else iFlag = 0;
#endif
	}
	va->push_back((*coordd)[1]);
	va->push_back((*coordd)[3]);
#if 0
	if(iFlag == 0)	{	ta->push_back(tStart[0]);	ta->push_back(tStart[1]);	}
	if(iFlag == 1)	{	ta->push_back(tEnd[0]);		ta->push_back(tEnd[1]);		}
#endif

	/* 计算 PrimitiveSet，不同方向要有不同的排列顺序。*/
	osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	unsigned int left_u = 0, left_d = 1;
	unsigned int right_u, right_d;
	for(unsigned int ii = 2; ii < va->size(); ii+=2)
	{
		right_u = ii; right_d = ii + 1;
		if(Direction == 1)
		{
			dui->push_back(left_u); dui->push_back(right_u); dui->push_back(left_d);
			dui->push_back(right_d); dui->push_back(left_d); dui->push_back(right_u);
		}
		if(Direction == 0)
		{
			dui->push_back(left_d); dui->push_back(right_u); dui->push_back(left_u); 
			dui->push_back(right_u); dui->push_back(left_d); dui->push_back(right_d); 
		}
		left_u = right_u; left_d = right_d;
	}
	gm2->addPrimitiveSet(dui.get());
	gd->removeDrawables(0, gd->getNumDrawables());
	gd->addDrawable(gm2.get());

#if 0
	osgDB::writeNodeFile(*gd, "F:\\g2.ive");
#endif

	/* 调整一下所有的顶点值 */
	/* 计算顶点高度的时候，要保证路面水平，也就是路面横向两个相对的路边点高度相同。 2012-10-18 */
	int i=0, iFlag2 = 0;
	double currElev;
	for(osg::Vec3Array::iterator v = va->begin(); v != va->end(); v++)
	{
		if(!iFlag2)
		{
			cart.x = v->x();  cart.y = v->y();  cart.z = v->z();
			CartToGeod(&cart, &geod);
			currElev = geod._elevation = GetElevation(geod._lon, geod._lat, 0.0);
			GeodToCart(&geod, &cart);
			v->set(cart.x, cart.y, cart.z);
			iFlag2 = 1;
		}else{
			cart.x = v->x();  cart.y = v->y();  cart.z = v->z();
			CartToGeod(&cart, &geod);
			geod._elevation = currElev;
			GeodToCart(&geod, &cart);
			v->set(cart.x, cart.y, cart.z);
			iFlag2 = 0;
		}
		
#if	0
	/* 暂不生成道路上的随机灯点。。2013-01-22 */

		/* 把这些顶点也加入到灯光序列中 */
		if(CurrentLampType == false)
		{
			if((hasLight) && (i%2))
			{
				SimPoint sp;
				sp.name = "USR_WHITE_LIGHTS_0_30000";
				sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
				sp.x    = cart.x;		sp.y	= cart.y;		sp.z	= cart.z;
				InsertCultPoint(sp);
			}
			i++;
		}
	/* 注释完毕 */
#endif	//0

	}
#if 0
	osgDB::writeNodeFile(*gd, "F:\\g3.ive");
#endif

	return gd.release();
}

/* 合并左边路和右边路，组合成一条大路 */
osg::ref_ptr<osg::Geode> MergeRoadFromLAndR(osg::ref_ptr<osg::Geode> gLRoad, osg::ref_ptr<osg::Geode> gRRoad, int Direction)
{
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	osg::ref_ptr<osg::Vec3Array> vRoad = new osg::Vec3Array;
#if 1
	osg::ref_ptr<osg::Vec2Array> tRoad = new osg::Vec2Array;
	osg::Vec2 tStart[2], tEnd[2];	tStart[0].set(0.0, 0.0);	tStart[1].set(1.0, 0.0);	tEnd[0].set(0.0, 1.0);		tEnd[1].set(1.0, 1.0);
#endif
	int iFlag = 0;
	
	/* 先取出左侧边路和右侧边路的顶点，交叉组合成一个顶点阵列 */
	osg::ref_ptr<osg::Geometry> lroad = dynamic_cast<osg::Geometry*> (gLRoad->getDrawable(0));
	osg::ref_ptr<osg::Geometry> rroad = dynamic_cast<osg::Geometry*> (gRRoad->getDrawable(0));
	osg::Vec3Array *vl = dynamic_cast<osg::Vec3Array *>(lroad->getVertexArray());
	osg::Vec3Array *vr = dynamic_cast<osg::Vec3Array *>(rroad->getVertexArray());
	
	/* 把左右顶点组各自里面的奇数点提出来，*/
	osg::Vec3Array::iterator iL = vl->begin();	iL++;
	vRoad->push_back(*iL); iL += 2;
	osg::Vec3Array::iterator iR = vr->begin();	iR++;
	vRoad->push_back(*iR); iR += 2;
#if 1
	tRoad->push_back(tStart[0]);	tRoad->push_back(tStart[1]);	iFlag++;
#endif
	
	int i=1;
	for(iL; iL != vl->end(); iL++)
	{
		if(i%2)
		{
			vRoad->push_back(*iL);
			vRoad->push_back(*iR);
#if 1
			if(iFlag == 0)	{	tRoad->push_back(tStart[0]);	tRoad->push_back(tStart[1]);	}
			if(iFlag == 1)	{	tRoad->push_back(tEnd[0]);		tRoad->push_back(tEnd[1]);		}
			if(iFlag == 0)	iFlag = 1; else iFlag = 0;
#endif
		}
		i++;	iR++;
	}
	
	/* 组合成一个紧凑的Geode，并增加纹理坐标，纹理坐标增加方法是每个四边形用一个整图。 */
	osg::ref_ptr<osg::Geometry>  gm = new osg::Geometry;
	gm->setVertexArray(vRoad);
#if 1
	gm->setTexCoordArray(0, tRoad);
#endif
	
	/* 计算 PrimitiveSet，不同方向要有不同的排列顺序。*/
	osg::ref_ptr<osg::DrawElementsUInt> dui = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES);
	unsigned int left_u = 0, left_d = 1;
	unsigned int right_u, right_d;
	for(unsigned int ii = 2; ii < vRoad->size(); ii+=2)
	{
		right_u = ii; right_d = ii + 1;
		if(Direction == 1)
		{
			dui->push_back(left_u); dui->push_back(right_u); dui->push_back(left_d);
			dui->push_back(right_d); dui->push_back(left_d); dui->push_back(right_u);
		}
		if(Direction == 0)
		{
			dui->push_back(left_d); dui->push_back(right_u); dui->push_back(left_u); 
			dui->push_back(right_u); dui->push_back(left_d); dui->push_back(right_d); 
		}
		left_u = right_u; left_d = right_d;
	}
	gm->addPrimitiveSet(dui.get());
	gd->removeDrawables(0, gd->getNumDrawables());
	gd->addDrawable(gm.get());
	
	return gd.release();
}


void CShpFile::GetCurrLinePointElevation(Geodesy *geod, int mode)
{
	if(mode == NEED_ELEVATION)
	{
		geod->_elevation = GetElevation(geod->_lon, geod->_lat, 0.0);
	}else
		geod->_elevation = 0.0;
}


osg::ref_ptr<osg::Group> CShpFile::DrawPLineElement(int mode, CMapLine* pPline, float roadWidth, int type)
{
	unsigned int test = 4;

	int i,j,k,iCount, isLastPointValid;
	CMapPoint *pPoint;
	CMapPoints *pPoints;
	CMapParts *pParts;
	Geodesy geod;
	Carton  cart;
	osg::Vec3 lastPoint;
	osg::Vec3 currPoint; 
	
	totallines.clear();
	MyPolyLine aline; 
	int currMode = NO_ELEVATION;
	
	if ( pPline == NULL )
		return NULL;
	
	/* 先将所有的线都提取出来，保存在totallines里面。*/
	for ( i = 0 ; i < pPline->GetCount() ; i++ )
	{
		aline.clear();
		pParts =  pPline->GetParts(i);
		isLastPointValid = 0;

		for ( j = 0 ; j < pParts->GetCount() ; j++)
		{
			pPoints = pParts->Get(j); 
			iCount = pPoints->GetCount();
			
			osg::Vec3 start;
			pPoint = pPoints->Get(0);

			geod._lon = pPoint->GetX();
			geod._lat = pPoint->GetY();
			SetMMLL(geod._lon, geod._lat);
			lastPoint.set(geod._lon, geod._lat, 0.0);
			isLastPointValid = USE_INSERT_VERTEX;
			
			GetCurrLinePointElevation(&geod, currMode);		//mode 
			GeodToCart(&geod, &cart);
			start.set(cart.x, cart.y, cart.z);
			
			for ( k = 1 ; k < iCount ; k++)
			{	
				pPoint = pPoints->Get(k); 
				osg::Vec3 end;

				geod._lon = pPoint->GetX();
				geod._lat = pPoint->GetY();
				SetMMLL(geod._lon, geod._lat);

				currPoint.set(geod._lon, geod._lat, 0.0);
				if(isLastPointValid)
				{
					/* 首先在上一点和当前点之间按照固定的距离插入若干中间点。*/
					std::vector<osg::Vec3> vec_ab;
					interpolation(lastPoint, currPoint, vec_ab, INSERTDISTANCE);
					vec_ab.push_back(currPoint);
					
					/* 把每个插入点都计算高程，并插入点队列中 */
					for(unsigned int m=0; m<vec_ab.size(); m++)
					{
						geod._lon = vec_ab[m].x(); geod._lat = vec_ab[m].y();
						GetCurrLinePointElevation(&geod, currMode);		//mode
						GeodToCart(&geod, &cart);
						end.set(cart.x, cart.y, cart.z);
	
						SimEdge se;
						se.A = start; se.B = end; se.x = k-1; se.y = k; se.z = 0;
						aline.push_back(se);
						start = end;
					}
				}else
				{
					GetCurrLinePointElevation(&geod, currMode);		//mode
					GeodToCart(&geod, &cart);
					end.set(cart.x, cart.y, cart.z);
	    			
					SimEdge se;
					se.A = start; se.B = end; se.x = k-1; se.y = k; se.z = 0;
					aline.push_back(se);
					start = end;
				}

				/* 插入点的起始点 */
				lastPoint = currPoint;
				isLastPointValid = USE_INSERT_VERTEX;
			}
		}// j
		totallines.push_back(aline);		
	}

	/* 根据线组创建道路及两边的树 */
	if(mode == NEED_ELEVATION)
	{
		osg::ref_ptr<osg::Group> gp;

		/* type == 0 : Road; type == 1 : Forest */
		if(type == 0)
			gp = CreateRoadsWithTree(roadWidth);
		if(type == 1)
			gp = CreateTreeLine(roadWidth);

		return gp.release();
	}
	else
		return NULL;
}

/************************************************************

  MapLayer

***********************************************************/

CMapLayer::CMapLayer()
{
	m_bVisible = TRUE;
	m_Valid = FALSE;
	m_csLayerName = _T("");
	m_pRender = new CMapRender;
    m_pRender->SetRenderType(SIMPLE_RENDER);
}

CMapLayer::~CMapLayer()
{
	if ( m_pRender != NULL )
		delete m_pRender;
}

BOOL CMapLayer::LoadData(CString& csFileName)
{
	CString csDbfName;
	if (m_shpFile.ReadShp(csFileName))
    {
		m_Valid = TRUE;
		return TRUE;
    }
	else
    {
		m_Valid = FALSE;
		return FALSE;
	}
}

CMapRectangle CMapLayer::GetExtent()
{
	return m_shpFile.GetExtent();
}

long CMapLayer::GetShpType()
{
	return m_shpFile.GetShpType();
}

BOOL CMapLayer::GetVisible()
{
	return m_bVisible;
}

void CMapLayer::SetVisible(bool bVisible)
{
	m_bVisible = bVisible;
}

void CMapLayer::SetLayerName(CString& csLayerName )
{
	m_csLayerName = csLayerName;
}

CString CMapLayer::GetLayerName()
{
	return m_csLayerName;
}

void CMapLayer::SetExtent(CMapRectangle& extent )
{
	m_shpFile.SetExtent(extent);
}

void CMapLayer::SetShpType(int lShpType )
{
	m_shpFile.SetShpType(lShpType);
}
int CMapLayer::SearchShape(CMapPoint& pt )
{
      if ( GetShpType() != POLYGON )
		  return -1;
	  else
      {
		 return m_shpFile.SearchShape(pt);
	  }
}
void CMapLayer::SetRender(CMapRender* pRender )
{
	if ( m_pRender != NULL )
		delete m_pRender;
	m_pRender = pRender;
}

osg::ref_ptr<osg::Group> CMapLayer::DrawLayer(int mode)
{
	osg::ref_ptr<osg::Group> gp;
	gp = m_shpFile.DrawShp(mode);
	return gp;
}

/************************************************************

 	MapLayers

***********************************************************/
CMapLayers::CMapLayers()
{
}

CMapLayers::~CMapLayers()
{
	Clear();
}

short CMapLayers::GetCount()
{
	return m_Layers.GetSize();
}

void CMapLayers::SetCount(short sCount)
{
}

CMapLayer* CMapLayers::GetAt(short sIndex)
{

	int iCount;
    CMapLayer *pLayer;

    iCount = m_Layers.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount)
		return NULL;
	pLayer = m_Layers.GetAt(sIndex);
	return pLayer;
}

BOOL CMapLayers::Add(CMapLayer* pMapLayer)
{
	if ( pMapLayer == NULL )
		return FALSE;
	if (!m_Layers.Add( pMapLayer ))
		return FALSE;
	return true;
}

void CMapLayers::Remove(short sIndex)
{
	int iCount;
    CMapLayer *pLayer;

    iCount = m_Layers.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount)
		return;

	pLayer = m_Layers.GetAt(sIndex);
	m_Layers.RemoveAt(sIndex);
	delete pLayer;
}

void CMapLayers::Clear()
{
   int i,iCount;
   CMapLayer *pLayer;

   iCount = m_Layers.GetSize() - 1;
   for ( i = iCount ; i >= 0   ; i-- )
   {
		pLayer = m_Layers.GetAt(i);
		delete pLayer;
   }
   m_Layers.RemoveAll();
}

void CMapLayers::MoveTo(short fromIndex, short toIndex)
{
	int iCount;
    CMapLayer *pLayer;

    iCount = m_Layers.GetSize() - 1;
	if ( fromIndex < 0 || fromIndex > iCount)
		return;
    if ( toIndex < 0 || toIndex > iCount)
		return;
	pLayer = m_Layers.GetAt(fromIndex);
	m_Layers.RemoveAt(fromIndex);
    m_Layers.InsertAt(toIndex, pLayer);
}

void CMapLayers::MoveToTop(short sIndex)
{
	int iCount;
    CMapLayer *pLayer;

    iCount = m_Layers.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount)
		return;
	pLayer = m_Layers.GetAt(sIndex);
	m_Layers.RemoveAt(sIndex);
    m_Layers.InsertAt(0, pLayer);
}

void CMapLayers::MoveToBottom(short sIndex)
{
	int iCount;
    CMapLayer *pLayer;

    iCount = m_Layers.GetSize() - 1;
	if ( sIndex < 0 || sIndex > iCount)
		return;
	pLayer = m_Layers.GetAt(sIndex);
	m_Layers.RemoveAt(sIndex);
    m_Layers.InsertAt(iCount, pLayer);
}

void CMapLayers:: GetAllExtent(CMapRectangle& rc)
{
	short i;
	CMapLayer *pLayer;

    for ( i = 0 ; i < m_Layers.GetSize() ; i++ )
	{
		pLayer = m_Layers.GetAt(i);
		CMapRectangle extent(pLayer->GetExtent());
		if ( i == 0 )
		{
			rc.SetLeft(extent.GetLeft());
			rc.SetRight(extent.GetRight());
			rc.SetTop(extent.GetTop());
			rc.SetBottom(extent.GetBottom());
		}
		else
        {
           if ( extent.GetLeft() < rc.GetLeft() )
			   rc.SetLeft(extent.GetLeft());
		   if ( extent.GetRight() > rc.GetRight() )
			   rc.SetRight(extent.GetRight());
           if ( extent.GetTop() < rc.GetTop() )
			   rc.SetTop(extent.GetTop());
		   if ( extent.GetBottom() > rc.GetBottom() )
			   rc.SetBottom(extent.GetBottom());
		}
     }
}

/************************************************************

 	Interface

***********************************************************/
void InitMMLL(CString fn)
{
	string filename = fn;
	currLayer.mShpName = filename;
	currLayer.minLong = currLayer.minLati = 180.0;
	currLayer.maxLong = currLayer.maxLati = -180.0;
}

void SetMMLL(double lon, double lat)
{
	if(lon > currLayer.maxLong) currLayer.maxLong = lon;
	if(lat > currLayer.maxLati) currLayer.maxLati = lat;
	if(lon < currLayer.minLong) currLayer.minLong = lon;
	if(lat < currLayer.minLati) currLayer.minLati = lat;
}

int LoadShpAndDrawIt(char *shpFn)
{
	int ret = 0;
	CString csFileName(shpFn);
	osg::ref_ptr<osg::Group> group;

	CultPoints.clear();			/* 点光源列表初始化。*/
	TreePoints.clear();			/* 点树列表初始化。*/
	CMapLayer* pLayer = new CMapLayer;
	pLayer->SetLayerName(csFileName); 
	if ( pLayer->LoadData(csFileName) > 0 )
	{
		InitMMLL(csFileName);	
		group = pLayer->DrawLayer(NEED_ELEVATION);
		if(CreateDataStruFromGroup(group))
			ret = 1;
		else
		{
			/* 如果当前生成树和灯，则返回2 */
			if((CultPoints.size() > 0)||(TreePoints.size() > 0))
				ret = 2;
		}
	}

	return ret;
}

int LoadShpAndCalculateMM(char *shpFn)
{
	int ret = 0;
	CString csFileName(shpFn);
	osg::ref_ptr<osg::Group> group;

	CultPoints.clear();			/* 点光源列表初始化。*/
	TreePoints.clear();			/* 点树列表初始化。*/
	CMapLayer* pLayer = new CMapLayer;
	pLayer->SetLayerName(csFileName); 
	if ( pLayer->LoadData(csFileName) > 0 )
	{
		InitMMLL(csFileName);
		group = pLayer->DrawLayer(NO_ELEVATION);
		ret = 1;
	}
	
	return ret;
}
