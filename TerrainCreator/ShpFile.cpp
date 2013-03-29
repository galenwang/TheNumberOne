#include "stdafx.h"
#include "d3dtypes.h"
#include "global.h"
#include "ShpFile.h"

#define USE_INSERT_VERTEX 1
#define INSERTDISTANCE 0.002

//�ⲿ�������ü�ȫ�ֱ���
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

//�����ö�����������
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
  ����:   �������
  ����: �������
  ����ֵ  ����

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
  ����: �ж�ָ�����Ƿ����߶���
  ����: ��p1 --- �߶���㡡p2 --- �߶��յ�
  ����ֵ  ���߶��� ����TRUE ���򷵻�FALSE

******************************************************************************/

bool   CMapPoint::IsPointInLine(CMapPoint& p1 , CMapPoint& p2 )
{
	double dblMulti;
	// �жϵ��Ƿ�����Ӿ��η�Χ��
	if ( m_dbX < min(p1.GetX() ,p2.GetX()) || m_dbX > max(p1.GetX() ,p2.GetX())
	 || m_dbY <  min(p1.GetY() ,p2.GetY()) || m_dbY > max(p1.GetY() ,p2.GetY()))
		return FALSE;
    //������
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
  ����:   �㵽�߾���
  ����: ����
  ����ֵ  ����

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
				//����㵽�߶���С����
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
  ����:    ����㵽�߶���С����
  ����: �� p1 --- �߶���㡡p2 --- �߶��յ�
  ����ֵ:  �ڵ㵽�߶���С����

******************************************************************************/
double CMapLine::ptToSegment(CMapPoint& pt,CMapPoint& ptStart,CMapPoint& ptEnd)
{
	double dblDist,dblX,dblY,k;
	CMapPoint pPlumb;

	if ( pt.IsPointInLine(ptStart, ptEnd) )
		return 0.0; //����ֱ����
    if ( fabs(ptEnd.GetX() - ptStart.GetX()) <= EP )
	{
		dblX = ptStart.GetX();
		dblY = pt.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.GetY() > min(ptStart.GetY(), ptEnd.GetY()) && pPlumb.GetY()
			< max(ptStart.GetY(), ptEnd.GetY()))
        {
			//�������߶η�Χ��
			return fabs(pPlumb.GetX() -pt.GetX());
        }
		else
        {   //�ж��߶ε��ĸ��˵��봹��ȽϽ�,Ȼ�������������
            if ( fabs(pPlumb.GetY() - ptStart.GetY()) < fabs(pPlumb.GetY() - ptEnd.GetY()))
			   dblDist = pt.Distance(ptStart);
			else
			   dblDist = pt.Distance(ptEnd);
        }
		//�߶�Ϊ��ֱ���
    } else if ( fabs(ptEnd.GetY() - ptStart.GetY()) <= EP )
    {
		//�߶�Ϊˮƽ���
		dblX = pt.GetX();
		dblY = ptStart.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.GetX() > min(ptStart.GetX(), ptEnd.GetX()) && pPlumb.GetX()
			< max(ptStart.GetX(), ptEnd.GetX()))
        {
			//�������߶η�Χ��
			return fabs(pPlumb.GetY() -pt.GetY());
        }
		else
        {   //�ж��߶ε��ĸ��˵��봹��ȽϽ�,Ȼ�������������
            if ( fabs(pPlumb.GetX() - ptStart.GetX()) < fabs(pPlumb.GetX() - ptEnd.GetX()))
			   dblDist = pt.Distance(ptStart);
			else
			   dblDist = pt.Distance(ptEnd);
        }
    }
	else
    {
		//�߶���б״̬
		k = (double)((ptEnd.GetY() - ptStart.GetY() ) /(ptEnd.GetX() - ptStart.GetX()));
		// ���㴹������
		dblX = (k*k*ptStart.GetX() + k*(pt.GetY()-ptStart.GetY())+pt.GetX())/(k*k+1);
		dblY = k*(dblX-ptStart.GetX()) + ptStart.GetY();
		pPlumb.SetX(dblX);
		pPlumb.SetY(dblY);
		if ( pPlumb.IsPointInLine(ptStart,ptEnd) )
        {
			//�������߶���
		    dblDist = pt.Distance(pPlumb);
        }
		else
        {
			//�ж��߶ε��ĸ��˵��봹��ȽϽ�,Ȼ�������������
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
  ����:         �жϵ��Ƿ��ڶ������
  ����:         �����
  �����        �ڷ���TRUE ���ڷ���FALSE
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
	//���ھ��η�Χ��ֱ�ӷ���
	if ( !m_Extent.IsPointIn(pt) )
		  return FALSE;
	//��һ��ͨ��pt�������
	dblTemp = pt.GetY();
	ptInfint.SetY( dblTemp );
	dblTemp = infinity;
	ptInfint.SetX(dblTemp);
	//��ø��϶���ε�ÿһ�������
	for ( i = 0 ; i < m_Polygon.GetSize() ; i++ )
    {
		pParts = (CMapParts*)m_Polygon.GetAt(i);
		//���һ������εĵ㼯��
		for ( j = 0 ; j < pParts->GetCount() ; j++ )
		{
			pPoints = (CMapPoints*)pParts->Get(j);
			//���ÿ���㼯�ϵĶ�������
			for ( k = 0 ; k < pPoints->GetCount() - 1 ; k++ )
			{
				ptFirst  = (CMapPoint*)pPoints->Get(k);
				ptSecond = (CMapPoint*)pPoints->Get(k+1);
				if (pt.IsPointInLine(*ptFirst,*ptSecond) )
					return TRUE; //�õ��ڶ���α���
                if ( ptSecond->GetY() == ptFirst->GetY() )
					continue;  //�Թ�ˮƽ��

				if ( ptFirst->IsPointInLine(  pt , ptInfint ))
				{
					//���������ཻ�ڱߵ�һ������
					if ( ptFirst->GetY() > ptSecond->GetY() )
						iNumber++;
				}
				else if( ptSecond->IsPointInLine( pt , ptInfint ))
				{
					//���������ཻ�ڱߵĵڶ�������
					if ( ptFirst->GetY() > ptSecond->GetY() )
					    iNumber++;

				} //�ж��Ƿ��й�������
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
  ����:  �ж�����ֱ���Ƿ��ཻ �ҽ��㲻�Ƕ���
  ����: p1 --- �߶���� p2 ---�߶��յ� ��p3 --- �߶���㡡p4 --- �߶��յ�
  ����ֵ  ��ֱ���� ����TRUE ���򷵻�FALSE
******************************************************************************/

BOOL CMapPolygon::isIntersect(CMapPoint& p1 , CMapPoint& p2 , CMapPoint& p3 , CMapPoint& p4 )
{
	double dblMulti,dblTmp1,dblTmp2;
	//�����ཻ
	if ( p1.IsEqual(p3) || p1.IsEqual(p4) || p2.IsEqual(p3) || p2.IsEqual(p4) )
		return FALSE;
    //�ж������߶���Ӿ����Ƿ��ཻ
	if ( min(p1.GetX(),p2.GetX()) > max(p3.GetX(),p4.GetX()) || max(p1.GetX(),p2.GetX())
		< min(p3.GetX(),p4.GetX()) || min(p1.GetY(),p2.GetY()) > max(p3.GetY(),p4.GetY())
		|| max(p1.GetY(),p2.GetY()) < min(p3.GetY(),p4.GetY()))
		return FALSE;
    //������
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
  ����:         ��DBF�ļ�
  ����:         �ļ���
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
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

	//�����ļ�
	if ( !fDbf.Open(csFileName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	m_bDbfOpen = TRUE;
	fDbf.Seek(0L, CFile::begin);
	//�����ļ�ͷ
	if ( fDbf.Read(&m_Header,sizeof(m_Header)) != sizeof(m_Header))
    	return FALSE;
	ulReocrdCount = m_Header.no_recs;
	ulLength = m_Header.head_len;
    ulRecLength = m_Header.rec_len;

	//�����ֶθ���
	sFieldCount = (ulLength - sizeof(DBF_HEADER)-1)/sizeof(FIELD_ELEMENT);
	iTemp = sFieldCount * sizeof(FIELD_ELEMENT) + 1;
	pszBuffer = new char[iTemp];
	if ( pszBuffer == NULL )
		return FALSE;
    //�����ֶ�������������(��ṹ)
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

		//�ж��ֶ�����
		if ( pField->cFieldType != 'N' && pField->cFieldType != 'F' )
        {
			pField->ucFieldLength += pField->ucFieldDecimal*256;
            pField->ucFieldDecimal = 0;
		}
		pOldField = pField;
       	m_TableDesc.Add(pField);

    }
    //�����¼����¼��������
    ReadRecord(0);

	delete []pszBuffer;
	return true;
}

/*************************************************
  ����:         �������ݼ�¼
  ����:         ��¼������(��0��ʼ)
  �����        ��
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
        return; //��Ч������ֵ

    if ( iCursorPos != lRecordID) //Ҫ��ȡ�ļ�¼δ�ڻ�����
    {
		//�����¼����ļ�ͷ��ƫ����
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
			//�Թ��ü�¼�Ƿ�ɾ������ֽ�pszBuffer+1
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
		Clear(); //��ջ����������¼�¼
		m_Fields.Add(pMapFields);
		delete []pszBuffer;
		iCursorPos = lRecordID;
	}
}

/*************************************************
  ����:         �ƶ�����¼��ͷ��
  ����:         ��
  �����        ��
*************************************************/
void CMapRecordSet::MoveFirst()
{
	 if ( !m_bDbfOpen ) //���ݿ��ļ�δ��
		 return;
	 bBOF = TRUE;
	 bEOF = FALSE;
	 iCursorPos = -1;
	 ReadRecord(0);
}

/*************************************************
  ����:         �ƶ�����¼��β��
  ����:         ��
  �����        ��
*************************************************/
void CMapRecordSet::MoveLast()
{
     if ( !m_bDbfOpen )
		 return;
	 bEOF = TRUE;
	 ReadRecord(m_Header.no_recs - 1);
}

/*************************************************
  ����:         �ƶ�����һ����¼
  ����:         ��
  �����        ��
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
  ����:         �ƶ�����һ����¼
  ����:         ��
  �����        ��
*************************************************/
void CMapRecordSet::MovePrev()
{
	if ( !m_bDbfOpen ) //���ݿ��ļ�δ��
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
  ����:         �ƶ�iNumRecords����¼
  ����:         �ƶ��ļ�¼�����ƶ����λ��
  �����        ��
*************************************************/

BOOL CMapRecordSet::Move(int iNumRecords , RECORDSTART Start )
{
	int iPos;
	/*if ( bEOF && iNumRecords > 0 ) //�Ѿ�����¼��ĩβ
		return FALSE;
	if ( bBOF && iNumRecords < 0 ) //�Ѿ�����¼��ͷ
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

	if ( iNumRecords > 0 ) //����ƶ�
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
  ����:         ��ռ�¼������
  ����:         ��
  �����        ��
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
  ����:         ����Shp��shx�ļ� ��ʽ ����ESRI������Ƥ��
  ����:         �ļ���(ȫ·��)
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
*************************************************/


BOOL CShpFile::ReadShp(CString& csFileName)
{
	int	  iTemp;
	CString csShxName;
	CFileException fe;
    SHPHEADER varHeader;

	//�����ļ�
	if ( !fShp.Open(csFileName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	bShpOpen = TRUE;

	//�������ļ�
	csShxName = csFileName.Left(csFileName.GetLength() - 3);

	csShxName = csShxName + "shx";
    if ( !fShx.Open(csShxName, CFile::modeRead|CFile::shareDenyWrite,&fe))
		return FALSE;
	bShxOpen = TRUE;

	TRY
    {
	    //�����ļ�ͷ ��100�ֽ�
		if ( fShp.Read(&varHeader , sizeof(SHPHEADER))!= sizeof(SHPHEADER))
			return FILE_READERR;

        iTemp = varHeader.iFileCode;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&iTemp);
		if ( iTemp != 9994 ) //�Ƿ���shp�ļ�
			return FILE_CODEERR;
        if ( varHeader.iVersion != FILE_VERSION ) //�ļ��汾�Ƿ���ȷ
			return FILE_VERSIONERR;

		//shp ����
		m_shpType = varHeader.iShpType;
		m_shpFileLength = varHeader.iFileLength;
	    if ( !m_bBigEndian )
			SwapWord(sizeof(int),&m_shpFileLength);

		//�������������η�Χ
		m_Extent.SetLeft(varHeader.dbXMin);
		m_Extent.SetRight(varHeader.dbXMax);
		m_Extent.SetTop(varHeader.dbYMin);
		m_Extent.SetBottom(varHeader.dbYMax);

		//�������ļ�ͷ ��100�ֽ�
		if ( fShx.Read(&varHeader , sizeof(SHPHEADER))!= sizeof(SHPHEADER))
			return FILE_READERR;
		iTemp = varHeader.iFileCode;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&iTemp);
		if ( iTemp != 9994 ) //�Ƿ���shx�ļ�
			return FILE_CODEERR;
        if ( varHeader.iVersion != FILE_VERSION ) //�ļ��汾�Ƿ���ȷ
			return FILE_VERSIONERR;
		m_shxFileLength = varHeader.iFileLength;
		if ( !m_bBigEndian )
			SwapWord(sizeof(int),&m_shxFileLength);
		//ͨ�������ļ��������ļ���¼���� �ļ�������ֵ��16λ��
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
  ����:         ����DBF�ļ���ʽ
  ����:         �ļ���(ȫ·��)
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
*************************************************/
BOOL CShpFile::ReadDBF(CString& csFileName)
{
	CString csDbfName;
	BOOL bResult;
	//������¼������
	pRecordSet = new CMapRecordSet;
	ASSERT ( pRecordSet != NULL );
	csDbfName = csFileName.Left(csFileName.GetLength() - 3);
	csDbfName = csDbfName + "dbf";
	//��DBF�ļ�
	bResult =  pRecordSet->openDBF(csDbfName);
    if ( !bResult )
		delete pRecordSet;
	return bResult;
}

/*************************************************
  ����:         ��õ�ǰ��ͼ�ļ������η�Χ
  ����:         ��
  �����        ��ͼ���η�Χ
*************************************************/

CMapRectangle CShpFile::GetExtent()
{
	return m_Extent;
}

/*************************************************
  ����:         ���õ�ͼ�ļ������η�Χ
  ����:         ��ͼ���η�Χ
  �����        ��
*************************************************/
void CShpFile::SetExtent(CMapRectangle& extent )
{
	m_Extent.SetLeft(extent.GetLeft());
    m_Extent.SetRight(extent.GetRight());
	m_Extent.SetTop(extent.GetTop());
    m_Extent.SetBottom(extent.GetBottom());
}

/*************************************************
  ����:         ����Shp��������
  ����:         shp ����
  �����        ��
*************************************************/
void CShpFile::SetShpType(int& iShpType)
{
	m_shpType = iShpType;
}

/*************************************************
  ����:         ����Shp��������
  ����:         ��
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
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
				//������¼
				for ( i = 1 ; i <= m_iRecordCount ; i++ )
				{
					iLength = SetRecordPos(i);
					if ( iLength <= 0 )
						return FALSE;

					//��ü�¼ͷ��Ϣ
					if (!GetRecordHeader(RecordHeader))
						return FALSE;
					//��¼���ݳ����Ƿ���shpʵ���Сһ��,������¼�����Ƿ����¼���ݳ���һ��
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
					//�����Ƿ�ƥ��
					if ( point.iShpType != m_shpType)
						return FALSE;
                    pPoint->SetX(point.dbX );
					pPoint->SetY(point.dbY );
					m_ObList.AddTail(pPoint);

				 }
			}
			break;
        case POLYLINE:

			pszBuffer = new char [MAX_BUFFER_SIZE]; //���仺����
			if ( pszBuffer == NULL )
				return FALSE;
			memset(pszBuffer , 0 , MAX_BUFFER_SIZE);

			//�����߼�¼
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
				//��ü�¼ͷ��Ϣ
				if (!GetRecordHeader(RecordHeader))
					return FALSE;
				if ( !GetShpInfo(shpIn))
					return FALSE;
				 if (shpIn.ishpType != POLYLINE )//���Ͳ�ƥ��
                 {
                      delete pLine;
					  continue;
				 }
				 if ( shpIn.iNumParts*sizeof(int) > MAX_BUFFER_SIZE )
                 {
					 //�����߶���������󻺳�������,���Ըö���
					 delete pLine;
					 continue;
                 }
				 //�����������ʵ�ʳ���
				 j = sizeof(SHPINFO) + shpIn.iNumParts*sizeof(int) ;
				 j += shpIn.iNumPoints*sizeof(SHPPOINT);

				 //�ж�ʵ�ʳ����Ƿ��������ļ��м�¼��һ��
				 if ( RecordHeader.iContentLength*2 != j )
                 {
				     delete pLine;
					 continue;
				 }
				 //����shp���η�Χ
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
				 //��������߶�����
				 if ( fShp.Read(pIParts,shpIn.iNumParts*4) != (uint)(shpIn.iNumParts*4))
				 {
					delete pLine;
					return FALSE;
				 }

				 //������洢��ռ�ֽ���
				 iLength = shpIn.iNumPoints*sizeof(SHPPOINT);
				 //��ʼ������������
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
					    k = pIParts[j];        //���ε�һ����������
					    m = shpIn.iNumPoints ; //��һ���ε�һ����������
                    }
              	    else
                    {
						k = pIParts[j];
						m = pIParts[j+1];
					}
					//�����i�εĶ���
				    for ( ; k < m ; k++)
                    {
						pPoint = new CMapPoint();
						if ( pPoint == NULL )
							return FALSE;

						//��Ҫ�������ݸ��»�����
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
						//��Ҫ�������ݸ��»�����
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
			pszBuffer = new char [MAX_BUFFER_SIZE]; //���仺����
			if ( pszBuffer == NULL )
				return FALSE;
			memset(pszBuffer , 0 , MAX_BUFFER_SIZE);

			//�������μ�¼
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
				//��ü�¼ͷ��Ϣ
				if (!GetRecordHeader(RecordHeader))
					return FALSE;
				if ( !GetShpInfo(shpIn))
					return FALSE;
				 if (shpIn.ishpType != POLYGON )//���Ͳ�ƥ��
                 {
                      delete pPolygon;
					  continue;
				 }
				 if ( shpIn.iNumParts*sizeof(int) > MAX_BUFFER_SIZE )
                 {
					 //���϶�����еĶ���θ���������󻺳�������,���Ըö���
					 delete pPolygon;
					 continue;
                 }
				 //�����������ʵ�ʳ���
				 j = sizeof(SHPINFO) + shpIn.iNumParts*sizeof(int) ;
				 j += shpIn.iNumPoints*sizeof(SHPPOINT);

				 //�ж�ʵ�ʳ����Ƿ��������ļ��м�¼��һ��
				 if ( RecordHeader.iContentLength*2 != j )
                 {
				     delete pPolygon;
					 continue;
				 }
				 //����shp���η�Χ
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
				 //���븴�϶���Ͷ�����
				 if ( fShp.Read(pIParts,shpIn.iNumParts*4) != (uint)(shpIn.iNumParts*4))
				 {
					delete pPolygon;
					return FALSE;
				 }

				 //������洢��ռ�ֽ���
				 iLength = shpIn.iNumPoints*sizeof(SHPPOINT);
				 //��ʼ������������
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
					    k = pIParts[j];        //���ε�һ����������
					    m = shpIn.iNumPoints ; //��һ���ε�һ����������
                    }
              	    else
                    {
						k = pIParts[j];
						m = pIParts[j+1];
					}
					//�����i�εĶ���
				    for ( ; k < m ; k++)
                    {
						pPoint = new CMapPoint();
						if ( pPoint == NULL )
							return FALSE;

						//��Ҫ�������ݸ��»�����
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
						//��Ҫ�������ݸ��»�����
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
  ����:         ����ÿ��shp��������ļ�ͷ��ƫ����
  ����:         ��¼����ֵ(���㿪ʼ)
  �����        ��shp�����������ļ��е�λ��
*************************************************/
int CShpFile::SetRecordPos( int iRecord )
{
    unsigned int iOffset,iTmp;
	SHXRECORD shxRD;

	if ( iRecord < 0 )
		return 0;
    //��������ļ���¼ƫ��������ļ�ͷ
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
  ����:         ���ÿ��shp�����¼��¼ͷ����Ϣ
  ����:         ��¼ͷ�ṹ����
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
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
  ����:         ���ÿ��shp����������Ϣ
  ����:         ������Ϣ�ṹ����
  �����        �ɹ�����TRUE ʧ�ܷ���FALSE
*************************************************/

BOOL CShpFile::GetShpInfo(SHPINFO& varInfo)
{

	if(fShp.Read(&varInfo,sizeof(varInfo))!= sizeof(varInfo))
		return FALSE;
	return TRUE;
}

/*************************************************
  ����:         ������������
  ����:         ���ݻ�����ָ�� ���������32K
                ���������Ҫ�ֶ�ζ�ȡ��Ҫ��ȡ�ĳ��ȡ�
				�Ƿ��Ѷ�ȡ���
  �����        ��ȡ���ݵ�ʵ�ʳ���
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
  ����:         ���ݵ�������shp�����Ƿ�ѡ��
  ����:         �����
  �����        ���ҵ�shp���������ֵ ����ֵ-��ʾδ���ҵ�
*************************************************/
int CShpFile::SearchShape(CMapPoint& pt )
{
	unsigned long iCount;
	POSITION pos;
    CMapPolygon *pPolygon;

	if ( GetShpType() != POLYGON ) //ֻ�ж϶���ζ���
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
  ����:         ����shp����
  ����:         
  �����        osg::Group ָ��
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
  ����:         ���Ƶ����
  ����:         �豸ָ�롢ͼ������ָ�롢����任�����ṹ����
  �����        ��
*************************************************/
std::string DeleteSpace(CString csValue)
{
	/* ���csValue��β�пո����������� */
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
	
	/* ���sp.name���ǲ��ǵƵ�û��Lod������ */
	GetWidthFromName(sp.name);
	if((nameParts[0] == "USR") && (nameParts[2] == "LIGHTS"))
	{
		if(nameParts.size() <= 4) sp.name = sp.name + "_0_30000";
		sp.name = sp.name + "_" + sp.tile;
		CultPoints.push_back(sp);
	}
	
	/* ����ǵ������ӵ�TreePoints���� */
	if((nameParts[0] == "Forest") || (nameParts[0] == "forest") || (nameParts[0] == "EvergreenForest"))
	{
		sp.name = getTreeName(sp.name);
#if 1
		/* �������򣬰�������һ��߶ȡ�����1.5�� */
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
		
		/* ȡ������ */
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

		/* ȡ��X */	/* ȡ��Y */
		DrawPointElement(pPoint, &sp);

		/* ���Ҽ�ֵ */
		SetMMLL(sp.lon, sp.lat);

		/* ȡ�ø߳� */
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
  ����:         ���Ƶ����ĵ����
  ����:         �豸ָ�롢��������ָ�롢
                ����任�����ṹ����
  �����        ��
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
  ����ϸ����������
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
  ����:         ���ݶ���ζ��󴴽�Group
  ����:         ����ζ���
  �����        osg::Group
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
		
		/* �����ҵ�FieldNameΪ ��NAME�� ����*/
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

		/* �жϲ�����ϸ���������� */
		std::string dtname = DeleteSpace(csValue);
		InsertDetailedTxt(dtname);

		/* �����ǰ������Airport����ô��ȡ�����ĸ߶� */
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

		/* �����ǰ������Lake����ô��ȡ����ĸ߶� */
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

		
		/* �����ǰPolygon�ǵ�������ȫ�ֱ��� */
		if ((csValue[0] == 'I') && (csValue[1] == 's') && (csValue[2] == 'l') && (csValue[3] == 'a') && (csValue[4] == 'n') && (csValue[5] == 'd')) 
			ifHasIslandInCurrentTerrain = 1;
		
		/* ����ǰPolygon�Ƿ��Ǻ� */
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
  ����:         ���ݶ���ζ��󴴽�Geode
  ����:         ����ζ���
  �����        osg::Geode
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
					/* ��������һ��͵�ǰ��֮�䰴�չ̶��ľ�����������м�㡣*/
					std::vector<osg::Vec3> vec_ab;
					interpolation(lastPoint, currPoint, vec_ab, INSERTDISTANCE);
					vec_ab.push_back(currPoint);
					
					/* ��ÿ������㶼����̣߳������������� */
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

				/* ��������ʼ�� */
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
  ����:         �����ߵ����ƣ����»��߻���������Ϊ�ֶα�־�������Ʒ�Ϊ���Σ�
  ����:         string name
  �����        int width
*************************************************/
float GetWidthFromName(std::string matname)
{
	char matpart[128], mat[256];
	int j = 0, iWidth;
	
	/* ��matname���»�����Ϊ�ָ������ֳ�����С�Σ�������nameParts���� */
	sprintf(mat, "%s", matname.c_str());
	nameParts.clear();
	for(unsigned int i = 0; i<strlen(mat); i++)
	{
		if((mat[i] != '_')&&(mat[i] != '-')) { matpart[j++] = mat[i];	}
		if((mat[i] == '_')||(mat[i] == '-')) { matpart[j] = 0;	nameParts.push_back(matpart);	j = 0; }
	}
	matpart[j] = 0;	nameParts.push_back(matpart);	j = 0;
	
	/* ����nameParts�ĵڶ�����Ϊ���ֵ */ 
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
  ����:         ���ݶ���ζ��󴴽�Group
  ����:         ����ζ���
  �����        osg::Group
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

		/* �����ҵ�FieldNameΪ ��NAME�� ����*/
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
		
		/* ����������֣����еõ�Ҫ����·�Ŀ��ֵ������0~100֮�䣬�涨����ΪRoad_???��???Ϊ���ֵ  */
		std::string dtname = DeleteSpace(csValue);
		float curWidth = GetWidthFromName(dtname);
		
		/* �����������жϣ���ǰ������·�������� */
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
  ����:         ����PolyLine���󴴽�Geode
  ����:         ����ζ���
  �����        osg::Geode
*************************************************/

/* ��֪һ���߶��ϵ�A, B���㣬����Ҫ��AB�߶ηֶεĶ������󣺿���A��������һ���ֶε㡣*/
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
		/* fabs(delx)��fabs(dely)��ȣ�����������0.0 */
		count = fabs(delx) / step + 1;
		delta = delx / count;

		iVec[0] = A[0] + delta;
		iVec[1] = (delta)*dely/delx + A[1];
		iVec[2] = (delta)*delz/delx + A[2];
	}
	
	return iVec;
}

/* ��ȡ�ռ�������A��B֮����ֱ��AB�ϵĵ�����(segs)�������ȵĵ㡣 */
void GetLightPots(osg::Vec3d A, osg::Vec3d B, float distance)
{
	double delx, dely, delz, step, delta;
	std::vector<osg::Vec3d> vec_list;
	int count;
	Geodesy geod;
	Carton  cart;
	
	/* ����A�� */
	SimPoint sp;
	sp.name = "USR_WHITE_LIGHTS_0_30000";
	cart.x = sp.x = A.x();	cart.y = sp.y = A.y();	cart.z = sp.z = A.z();
	CartToGeod(&cart, &geod);
	sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
	InsertCultPoint(sp);
	
	/* AB֮�����segs�ݣ�����㶼���뵽vec_list���� */
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
		/* fabs(delx)��fabs(dely)��ȣ�����������0.0 */
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
	
	/* ��vec_list��������е㶼����CultPoint */
	for(unsigned int i=0; i<vec_list.size(); i++)
	{
		osg::Vec3d C = vec_list[i];
		
		sp.name = "USR_WHITE_LIGHTS_0_30000";
		cart.x = sp.x = C.x();	cart.y = sp.y = C.y();	cart.z = sp.z = C.z();
		CartToGeod(&cart, &geod);
		sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
		InsertCultPoint(sp);
	}
	
	/* ������B�� */
	sp.name = "USR_WHITE_LIGHTS_0_30000";
	cart.x = sp.x = B.x();	cart.y = sp.y = B.y();	cart.z = sp.z = B.z();
	CartToGeod(&cart, &geod);
	sp.lon  = geod._lon;	sp.lat  = geod._lat;	sp.elev = geod._elevation;
	InsertCultPoint(sp);	
}

/* �ڵ�ǰ���ϲ���Ƶ� */
double totalOffset = 0.0;
void InsertLightPotsOnCurrEdge(std::vector<GeoEdge>::iterator v, float distance)
{
	double currDist  = v->dist;							//��ǰ�����ߵĳ��ȡ�
	double validDist;
	int	   numLights;
	osg::Vec3d A = v->Ac, B = v->Bc;

	if(totalOffset > distance / 10.0)
	{
		/* һ�������v�߶��Ͽ���A�ˣ���B�˾���Ϊ(distance - totalOffset)�ĵ��ֵ: midA */
		osg::Vec3d midA = GetOffsetVertex(A, B, (distance - totalOffset));
			
		/* ���������midA��B֮���ܲ�����ٵƵ� */
		validDist = currDist - (distance - totalOffset);
		numLights = (int)( validDist / distance);
		
		/* �����������numLights��֮��ʣ�¶��پ��� */
		double currOffset = validDist - numLights * distance;
		if(currOffset > distance / 10.0)
		{
			/* �ȼ����v�߶��Ͽ���B�ˣ���A�˾���ΪcurrOffset�ĵ��ֵ: midC */
			osg::Vec3d midB = GetOffsetVertex(B, A, currOffset);
				
			/* �����߶�A<->midC֮�����numLights�ݣ�����·�Ƶ��ˣ���ת�浽CultPoints*/
			GetLightPots(midA, midB, distance);
			
			/* ʣ�µ�һ��ת��ΪtotalOffset������һ���߶�ʹ�� */
			totalOffset = currOffset;
		}else{
			/* currOffset ���Բ��� */
			/* ֱ�����߶�A<->B֮�����numLights�ݣ�����·�Ƶ��ˣ���ת�浽CultPoints*/
			GetLightPots(midA, B, distance);
		}
	}else{
		/* ���Բ��� totalOffset */
		/* һ�������A��B֮���ܲ�����ٵƵ� */
		validDist = currDist;
		numLights = (int)( validDist / distance);
		
		/* �����������numLights��֮��ʣ�¶��پ��� */
		double currOffset = validDist - numLights * distance;
		if(currOffset > distance / 10.0)
		{
			/* �ȼ����v�߶��Ͽ���B�ˣ���A�˾���ΪcurrOffset�ĵ��ֵ: midC */
			int  numOffsetSeg = (int)(validDist / currOffset);
			osg::Vec3d midB = GetOffsetVertex(B, A, numOffsetSeg);
				
			/* �����߶�A<->midC֮�����numLights�ݣ�����·�Ƶ��ˣ���ת�浽CultPoints*/
			GetLightPots(A, midB, distance);
			
			/* ʣ�µ�һ��ת��ΪtotalOffset������һ���߶�ʹ�� */
			totalOffset = currOffset;
		}else{
			/* currOffset ���Բ��� */
			/* ֱ�����߶�A<->B֮�����numLights�ݣ�����·�Ƶ��ˣ���ת�浽CultPoints*/
			GetLightPots(A, B, distance);
		}
	}
}

/* �ڱ�������һ���ļ�����ɵƵ� */
int CreateLightPotOnSideLine(float distance)
{
	totalOffset = 0.0;
	/* ���ݾ�γ�߸�ʽ�ı߼������� */
	std::vector<GeoEdge>::iterator v = SideGeoLine.begin();
	while(v != SideGeoLine.end())
	{
		double currDist  = v->dist;							//��ǰ�����ߵĳ��ȡ�
		int	   numLights = (int)( (currDist + totalOffset) / distance);		//��ǰ������Ӧ���еĵƵ���
		if(numLights > 0)
		{
			InsertLightPotsOnCurrEdge(v, distance);
		}else{
			totalOffset += currDist;
		}
		
		/* ��һ���� */
		v++;
	}
	
	return 1;	
}

/* ����һ���̻��� */
osg::ref_ptr<osg::Group> CreateTreeLine(float width)
{
	/* �������е�����ѭ�� */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	for(std::vector<MyPolyLine>::iterator tt = totallines.begin(); tt != totallines.end(); tt++)
	{
		/* ��� */
		/* ����֪���������һ���̻��� */
		MyPolyLine mpl = *tt;
		osg::ref_ptr<osg::Geode> gLTree = CreateRoadWithPolyLine(&mpl, width/2.0, 0, false);
		if(gLTree.valid())
		{
			gLTree->setName(DEFINE_FOREST_STRING);
			gp->addChild(gLTree.get());
		}

		/* �ұ� */
		/* ����֪���ұ���һ���̻��� */
		osg::ref_ptr<osg::Geode> gRTree = CreateRoadWithPolyLine(&mpl, width/2.0, 1, false);
		if(gRTree.valid())
		{
			gRTree->setName(DEFINE_FOREST_STRING);
			gp->addChild(gRTree.get());
		}
	}
	
	return gp.release();
}

/* ����һ�����ߴ��̻�����������·�Ƶ�·���Ļ���Ϣ */
osg::ref_ptr<osg::Group> CreateRoadsWithTree(float width)
{
	/* �������е�����ѭ�� */
	osg::ref_ptr<osg::Group> gp = new osg::Group;
	for(std::vector<MyPolyLine>::iterator tt = totallines.begin(); tt != totallines.end(); tt++)
	{
		/* ��� */
		/* �Ȱ���֪������һ��· */
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
		
		/* �������Geode�ҵ�ƽ���Ǹ��ߣ���·����ıߣ�������SideLine���� */
		if(gLRoad.valid())	CreateSideLineWithCurrentRoad(gLRoad);
		
#if 0
	/* �ݲ����ɵ�·�ϵ�����Ƶ㡣��2013-01-22 */
		/* ����·����ı���SideLine��������ȼ���ĵƵ� */
		if((mLightPotDistance > 0.0) && (CurrentLampType == true))
		{
			if(CreateLightPotOnSideLine(mLightPotDistance))
			{
			}
		}
	/* ע����� */	
#endif	//0
					
#if 0
	/* �ݲ����ɵ�·�ߵ��̻�������2013-01-23 */
		/* ���ݶԱ�����������̻��� */
		osg::ref_ptr<osg::Geode> gLeft  = CreateRoadWithPolyLine(&SideLine, 20.0, 0, false);
		if(gLeft.valid())
		{
			gLeft->setName(DEFINE_FOREST_STRING);
			gp->addChild(gLeft.get());
		}
	/* ע����� */	
#endif	//0

		/* �ұ� */
		/* �Ȱ���֪������һ��· */
		osg::ref_ptr<osg::Geode> gRRoad = CreateRoadWithPolyLine(&mpl, width/2.0, 1, true);
		if(gRRoad.valid())
		{
			//osgDB::writeNodeFile(*gLRoad, "F:\\grroad.osg");
			//gRRoad->setName(DEFINE_ROAD_STRING);
			//gp->addChild(gRRoad.get());
		}
		
		/* �������Geode�ҵ�ƽ���Ǹ��ߣ���·����ıߣ�������SideLine���� */
		if(gRRoad.valid())	CreateSideLineWithCurrentRoad(gRRoad);

#if 0
	/* �ݲ����ɵ�·�ϵ�����Ƶ㡣��2013-01-22 */
		/* ����·����ı���SideLine��������ȼ���ĵƵ� */
		if((mLightPotDistance > 0.0) && (CurrentLampType == true))
		{
			if(CreateLightPotOnSideLine(mLightPotDistance))
			{
			}
		}
	/* ע����� */	
#endif	//0
				
#if 0
	/* �ݲ����ɵ�·�ߵ��̻�������2013-01-23 */
		/* ���ݶԱ��������ұ��̻��� */
		osg::ref_ptr<osg::Geode> gRight  = CreateRoadWithPolyLine(&SideLine, 20.0, 1, false);
		if(gRight.valid())
		{
			gRight->setName(DEFINE_FOREST_STRING);
			gp->addChild(gRight.get());
		}
	/* ע����� */	
#endif	//0
		
		/* �ϲ����·���ұ�·����ϳ�һ����· */
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

/* ���ݵ�ǰ��·�����ߣ�����·�ı��� */
void CreateSideLineWithCurrentRoad(osg::ref_ptr<osg::Geode> gd)
{
	Geodesy geodA, geodB, gSL;
	Carton  cartA, cartB, cSL;
	osg::Vec3 st, ed;

	SideLine.clear();	SideGeoLine.clear();
	osg::Geometry *g = dynamic_cast<osg::Geometry *> (gd->getDrawable(0));
	osg::Vec3Array *va = dynamic_cast<osg::Vec3Array *>(g->getVertexArray());
	
	/* �������������������Ū��һ��MyPolyLine */
	osg::Vec3Array::iterator v = va->begin();	v++;
	osg::Vec3 start; start.set(*v);	 			v+=2;
	osg::Vec3 startGeo; startGeo.set(start);

	/*��ǰ��תΪ���θ߶� */
	cSL.x = start.x();	cSL.y = start.y(); cSL.z = start.z();	CartToGeod(&cSL, &gSL);	gSL._elevation = 0.0;
	GeodToCart(&gSL, &cSL);	start.x() = cSL.x;	start.y() = cSL.y;	start.z() = cSL.z;

	int i=1;
	for(v; v != va->end(); v++)
	{
		if(i%2)
		{
			/* ����ֱ������ϵ�ı� */
			osg::Vec3 end; 	end.set(*v);
			osg::Vec3 endGeo; endGeo.set(end);

			/*��ǰ��תΪ���θ߶� */
			cSL.x = end.x();	cSL.y = end.y(); cSL.z = end.z();	CartToGeod(&cSL, &gSL);	gSL._elevation = 0.0;
			GeodToCart(&gSL, &cSL);	end.x() = cSL.x;	end.y() = cSL.y;	end.z() = cSL.z;

			SimEdge se;
			se.A = start; se.B = end; se.x = 2*i - 1; se.y = 2*i + 1; se.z = 0;
			SideLine.push_back(se);	
			
			/* ���㾭γ������ϵ�ı� */
			cartA.x = startGeo.x();  	cartA.y = startGeo.y();  	cartA.z = startGeo.z();
			CartToGeod(&cartA, &geodA);
			cartB.x = endGeo.x();  	cartB.y = endGeo.y();  	cartB.z = endGeo.z();
			CartToGeod(&cartB, &geodB);
			GeoEdge ge; ge.A = geodA; ge.B = geodB; ge.Ac = start; ge.Bc = end;	ge.distance();
			SideGeoLine.push_back(ge);

			/* ��һ���� */
			start = end;	startGeo = endGeo;
		}
		i++;
	}

	return;
}

/* ���ݵ�ǰ��·�����ߣ�����һ��·�Ķ�������У���·�棬������·��������·�桱 */
/* ����� line��Ӧ���Ǻ��θ߶�Ϊ0.0���ߡ������ص�Geode, ���Ǽ����˺��θ߶ȵ��ߡ� */
osg::ref_ptr<osg::Geode> CreateRoadWithPolyLine(MyPolyLine *line, float width, int Direction, bool hasLight)
{
	Geodesy geod;
	Carton  cart;

	if(line->size() == 0) return NULL;

	/* ������֪��һ���ߣ������һ�������ľ��Σ����������ͬһ��ı���������֪�ߣ���������εĿ����ȡ�*/
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


	/* ��������������֮��Ĺ����� */
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

	/* ��ϳ�һ�����յ�Geode���������������꣬�����������ӷ�����ÿ���ı�����һ����ͼ�� */
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

	/* ���� PrimitiveSet����ͬ����Ҫ�в�ͬ������˳��*/
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

	/* ����һ�����еĶ���ֵ */
	/* ���㶥��߶ȵ�ʱ��Ҫ��֤·��ˮƽ��Ҳ����·�����������Ե�·�ߵ�߶���ͬ�� 2012-10-18 */
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
	/* �ݲ����ɵ�·�ϵ�����Ƶ㡣��2013-01-22 */

		/* ����Щ����Ҳ���뵽�ƹ������� */
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
	/* ע����� */
#endif	//0

	}
#if 0
	osgDB::writeNodeFile(*gd, "F:\\g3.ive");
#endif

	return gd.release();
}

/* �ϲ����·���ұ�·����ϳ�һ����· */
osg::ref_ptr<osg::Geode> MergeRoadFromLAndR(osg::ref_ptr<osg::Geode> gLRoad, osg::ref_ptr<osg::Geode> gRRoad, int Direction)
{
	osg::ref_ptr<osg::Geode> gd = new osg::Geode;
	osg::ref_ptr<osg::Vec3Array> vRoad = new osg::Vec3Array;
#if 1
	osg::ref_ptr<osg::Vec2Array> tRoad = new osg::Vec2Array;
	osg::Vec2 tStart[2], tEnd[2];	tStart[0].set(0.0, 0.0);	tStart[1].set(1.0, 0.0);	tEnd[0].set(0.0, 1.0);		tEnd[1].set(1.0, 1.0);
#endif
	int iFlag = 0;
	
	/* ��ȡ������·���Ҳ��·�Ķ��㣬������ϳ�һ���������� */
	osg::ref_ptr<osg::Geometry> lroad = dynamic_cast<osg::Geometry*> (gLRoad->getDrawable(0));
	osg::ref_ptr<osg::Geometry> rroad = dynamic_cast<osg::Geometry*> (gRRoad->getDrawable(0));
	osg::Vec3Array *vl = dynamic_cast<osg::Vec3Array *>(lroad->getVertexArray());
	osg::Vec3Array *vr = dynamic_cast<osg::Vec3Array *>(rroad->getVertexArray());
	
	/* �����Ҷ��������������������������*/
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
	
	/* ��ϳ�һ�����յ�Geode���������������꣬�����������ӷ�����ÿ���ı�����һ����ͼ�� */
	osg::ref_ptr<osg::Geometry>  gm = new osg::Geometry;
	gm->setVertexArray(vRoad);
#if 1
	gm->setTexCoordArray(0, tRoad);
#endif
	
	/* ���� PrimitiveSet����ͬ����Ҫ�в�ͬ������˳��*/
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
	
	/* �Ƚ����е��߶���ȡ������������totallines���档*/
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
					/* ��������һ��͵�ǰ��֮�䰴�չ̶��ľ�����������м�㡣*/
					std::vector<osg::Vec3> vec_ab;
					interpolation(lastPoint, currPoint, vec_ab, INSERTDISTANCE);
					vec_ab.push_back(currPoint);
					
					/* ��ÿ������㶼����̣߳������������� */
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

				/* ��������ʼ�� */
				lastPoint = currPoint;
				isLastPointValid = USE_INSERT_VERTEX;
			}
		}// j
		totallines.push_back(aline);		
	}

	/* �������鴴����·�����ߵ��� */
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

	CultPoints.clear();			/* ���Դ�б��ʼ����*/
	TreePoints.clear();			/* �����б��ʼ����*/
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
			/* �����ǰ�������͵ƣ��򷵻�2 */
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

	CultPoints.clear();			/* ���Դ�б��ʼ����*/
	TreePoints.clear();			/* �����б��ʼ����*/
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
