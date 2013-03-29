#include "stdafx.h"
#include "bucket.h"

/************************************************************
   
   MapPoint
   
***********************************************************/
#ifndef _MAP_POINT_H_
#define _MAP_POINT_H_

class CMapPoint:public CObject {

    DECLARE_DYNAMIC(CMapPoint)
public:
	CMapPoint();
    CMapPoint(CMapPoint& pt);  
	virtual ~CMapPoint();
public:
	void   SetX(double& dbX ) { m_dbX = dbX;};
	void   SetY(double& dbY ) { m_dbY = dbY;};
	void   SetIndex(unsigned int uiIndex ) { m_uiIndex = uiIndex ; };
	unsigned int GetIndex() { return m_uiIndex; }; 
	void   SetStatus(bool& bStatus) { m_bStatus = bStatus;};
	double GetX() { return m_dbX;};
	double GetY() { return m_dbY;};
	double Distance(CMapPoint& pt );
  	bool   GetStatus() { return m_bStatus;};
	bool   IsEqual(CMapPoint& pt );
	bool   IsPointInLine(CMapPoint& p1 , CMapPoint& p2 );
    
private:
	bool   m_bStatus;    //0---非选中状态, 1---选中状态
        
	unsigned int m_uiIndex;  //索引值
	double m_dbX;
	double m_dbY;
	CString m_csID;
};

#endif //_MAP_POINT_H_

/************************************************************
  
  MapRectangle
  
***********************************************************/

#ifndef _MAP_RECTANGLE_H_
#define _MAP_RECTANGLE_H_

class CMapRectangle {

public:
	CMapRectangle();
	CMapRectangle(CMapRectangle& MapRectangle );
	~CMapRectangle();
//attribute
public:
	void   SetLeft(double dbLeft) { m_dbLeft = dbLeft;};
	void   SetRight(double dbRight) { m_dbRight = dbRight;};
	void   SetTop(double dbTop ) { m_dbTop = dbTop;};  
	void   SetBottom(double dbBottom) { m_dbBottom = dbBottom; }; 
	double GetLeft()   { return m_dbLeft; };
	double GetRight()  { return m_dbRight; }
	double GetTop()    { return m_dbTop; };
	double GetBottom() { return m_dbBottom; };
	double GetWidth()  { return fabs(m_dbRight-m_dbLeft);};
	double GetHeigh()  { return fabs(m_dbBottom-m_dbTop);};
//operations
public:
	BOOL IsPointIn(CMapPoint& Point);
	BOOL IsInsercet(CMapRectangle& rc);
private:
	double m_dbLeft;
	double m_dbRight;
	double m_dbTop;
	double m_dbBottom;
};
#endif //_MAP_RECTANGLE_H_


/************************************************************
  
  MapPoints
  
***********************************************************/
#ifndef _MAP_POINTS_H_
#define _MAP_POINTS_H_

class CMapPoints:public CObject {

	DECLARE_DYNAMIC(CMapPoints)
public:
	CMapPoints();
    CMapPoints(CMapPoints& points);
	~CMapPoints();

//Attribute
public:
    long GetCount();
	CMapRectangle GetExtent();
	void SetExtent(CMapRectangle& exent);
//operations
public:
	CMapPoint* Get(long lIndex);
	void Add(CMapPoint* pPoint);
    void Set(long lIndex , CMapPoint* pPoint);
	void Remove(long lIndex);
	void Insert(long lndex , CMapPoint* pPoint);
	void Clear();
private:
	CMapRectangle m_Rectangle;
	CArray<CMapPoint*,CMapPoint*> m_Points; //存储点对象
};

#endif //_MAP_POINTS_H_


/************************************************************
  
  MapParts
  
***********************************************************/

#ifndef _MAP_PARTS_H_
#define _MAP_PARTS_H_

class CMapParts:public CObject {
	
	DECLARE_DYNAMIC(CMapParts)
public:
	CMapParts();
    CMapParts(CMapParts& Parts);
	~CMapParts();

//attribute
public:
	long GetCount();

//operation
public:
	void Add(CMapPoints* pPoints);
	void Set(long lindex, CMapPoints* pPoints);
	void Remove(long lindex);
	void Insert(long lindex, CMapPoints* pPoints);
	CMapPoints* Get(long lindex);
	void Clear();

private:
	CArray<CMapPoints*,CMapPoints*> m_Parts; //存储点集合对象

};

#endif //_MAP_PARTS_H_

/************************************************************
  
  MapLine
  
***********************************************************/
#ifndef _MAP_LINE_H_
#define _MAP_LINE_H_

class CMapLine:public CObject 
{
    
	DECLARE_DYNAMIC(CMapLine)
public:	
	CMapLine();
	CMapLine(CMapLine& mapline );
	~CMapLine();

//Attributes
public:
   long GetCount();
   CMapRectangle GetExtent();
   void SetExtent(CMapRectangle& exent);
   CMapParts* GetParts(long lIndex);
   double GetLength();
   double Distance(CMapPoint& pt );
   double ptToSegment(CMapPoint& pt,CMapPoint& ptStart,CMapPoint& ptEnd);

//operation
public:
	void Add(CMapParts* pParts);
	void Set(long lIndex , CMapParts* pParts);
	void SetIndex(unsigned int uiIndex ) { m_uiIndex = uiIndex ; }; 
	unsigned int GetIndex() { return m_uiIndex; }; 
	void Remove(long lIndex);
	void Insert(long lndex , CMapParts* pParts);
    void Clear();

private:
	unsigned int m_uiIndex;  //索引值
	CString m_csID;
	CMapRectangle m_Extent;
	CArray<CMapParts*,CMapParts*> m_Line; //存储点集合对象

};

#endif //_MAP_LINE_H_

/************************************************************
  
  MapPolygon
  
***********************************************************/

#ifndef _MAP_POLYGON_H_
#define _MAP_POLYGON_H_

class CMapPolygon:public CObject {

	DECLARE_DYNAMIC(CMapPolygon)
public:	
	CMapPolygon();
	CMapPolygon(CMapPolygon& mappolygon );
	~CMapPolygon();

//Attributes
public:
   long GetCount();
   CMapRectangle GetExtent();
   void SetExtent(CMapRectangle& exent);
   CMapParts* GetParts(long lIndex);
   double GetArea();

//operation
public:
	void Add(CMapParts* pParts);
	void Set(long lIndex , CMapParts* pParts);
	void SetIndex(unsigned int uiIndex ) { m_uiIndex = uiIndex ; }; 
	unsigned int GetIndex() { return m_uiIndex; }; 
	void Remove(long lIndex);
	void Insert(long lndex , CMapParts* pParts);
    void Clear();
	BOOL IsPointIn(CMapPoint& pt );
protected:
	BOOL isIntersect(CMapPoint& p1 , CMapPoint& p2 , CMapPoint& p3 , CMapPoint& p4 );
private:
	unsigned int m_uiIndex;  //索引值
	CMapRectangle m_Extent;
	CArray<CMapParts*,CMapParts*> m_Polygon; //存储点集合对象

};

#endif //_MAP_POLYGON_H_

/************************************************************
  
  DbfFile
  
***********************************************************/
#ifndef _DBFFILE_H_
#define _DBFFILE_H_

/* DBF文件头结构 */ 
typedef struct dbf_header 
{ 
    char vers;                        /* 版本标志*/ 
    unsigned char  yy,mm,dd;          /* 最后更新年、月、日 */ 
    unsigned long  no_recs;           /* 文件包含的总记录数 */ 
    unsigned short head_len;          /*文件头长度*/
	unsigned short rec_len;            /*记录长度 */ 
    char reserved[20];                 /* 保留 */ 
} DBF_HEADER; 

/*字段描述结构*/
typedef struct field_element {  
 
	char szFieldName[11];   /* 字段名称 */ 
    char cFieldType;        /* 字段类型 */ 
    unsigned long ulOffset; /* 偏移量 */ 
	unsigned char ucFieldLength; /* 字段长度 */ 
    unsigned char ucFieldDecimal; /* 浮点数整数部分长度 */ 
    char reserved1[2]; /* 保留 */ 
    char dbaseiv_id; /* dBASE IV work area id */ 
    char reserved2[10]; 
    char cProductionIndex;  
} FIELD_ELEMENT; 

#endif //_DBFFILE_H_

/************************************************************
  
  MapField
  
***********************************************************/

#ifndef _MAP_FIELD_H_
#define _MAP_FIELD_H_

typedef enum 
{
	
	fdString,    //字符串类型
	fdInteger,   //整型
	fdDouble,    //浮点型
	fdInvaild    //未知类型

} DBFFIELDTYPE;

class CMapField {

public:
	CMapField();
	CMapField(CMapField& field);
	~CMapField();

//Attribute
public:
	CString GetName();
	void SetName(LPCTSTR);
	long GetType();
	void SetType(long);
	CString GetValueAsString();
	void SetValueAsString(LPCTSTR);
	VARIANT GetValue();
	void SetValue(const VARIANT&);
	//VARIANT Get_Value();
	//void Set_Value(const VARIANT&);
private:
	CString m_csFieldName;
	long    m_lFieldType;
    VARIANT m_varValue;
	CString m_csValue;

};
#endif //_MAP_FIELD_H_


/************************************************************
  
  MapFields
  
***********************************************************/


#ifndef _MAP_FIELDS_H_
#define _MAP_FIELDS_H_

class CMapFields {

public:
	CMapFields();
    CMapFields(CMapFields& fields ); 
    ~CMapFields();

//Attribute
public:
	short GetCout();

//operation
public:
	void Add(CMapField* pField);
	void Remove(short sindex);
	void Insert(short sindex, CMapField* pField);
	CMapField *GetField(short sIndex);
	void Clear();

private:
	CArray<CMapField*,CMapField*> m_Fields; //存储字段
};

#endif //_MAP_FIELDS_H_


/************************************************************
  
  MapTableDesc
  
***********************************************************/
#ifndef _MAP_TABLEDESC_H_
#define _MAP_TABLEDESC_H_
#define INVALID_VALUE -1

class CMapTableDesc {

public:
	CMapTableDesc();
    CMapTableDesc(CMapTableDesc& tblDesc);
	~CMapTableDesc();

public:
	short GetFieldCount();
	void SetFieldCount(short);
	
// Operations
public:
	CString GetFieldName(short index);
	void SetFieldName(short index, LPCTSTR lpszNewValue);
	long GetFieldType(short index);
	void SetFieldType(short index, long nNewValue);
	short GetFieldPrecision(short index);
	void SetFieldPrecision(short index, short nNewValue);
	short GetFieldLength(short index);
	void SetFieldLength(short index, short nNewValue);
	short GetFieldScale(short index);
	void SetFieldScale(short index, short nNewValue);
    FIELD_ELEMENT *GetDesc(short sIndex);
	void Add(FIELD_ELEMENT* pField);
	void Remove(short sindex);
	void Insert(short sindex, FIELD_ELEMENT* pField);
	void Clear();

private:
	CArray<FIELD_ELEMENT*,FIELD_ELEMENT*> m_fieldsDesc; //存储表结构

};
#endif //_MAP_TABLEDESC_H_


/************************************************************
  
  MapRecordSet
  
***********************************************************/

#ifndef _MAP_RECORDSET_H_
#define _MAP_RECORDSET_H_
#define MAX_CACH_SIZE 100   //最大缓存记录数
class CMapFields;

//记录集移动的位置
typedef enum {
	BookmarkCurrent,
	BookmarkFirst,
	BookmarkLast

} RECORDSTART;
class CMapRecordSet {

public:
    CMapRecordSet();
	~CMapRecordSet();

//ATTRIBUTE
public:
	long GetRecordCount();
	CMapFields* GetFields(long sIndex);
	CMapTableDesc* GetTableDesc();
	BOOL GetBOF();
	BOOL GetEOF();
	int  GetCacheSize();
    BOOL SetCacheSize(int& CacheSize);   
//opeeration
public:
	BOOL openDBF(CString& csFileName);
	void MoveFirst();
	void MoveLast();
	void MoveNext();
	void MovePrev();
	BOOL Move(int iNumRecords , RECORDSTART Start );
private:
	void Clear();
	void ReadRecord(unsigned long lRecordID);
private:
	CArray<CMapFields*,CMapFields*> m_Fields; //记录集缓冲区
	CMapTableDesc  m_TableDesc;
	int m_CacheSize;
	unsigned long iCursorPos;                //游标当前位置
	BOOL bBOF,bEOF;

private:
	 DBF_HEADER m_Header;          //存储DBF文件头
	 BOOL  m_bDbfOpen;             //数据库文件是否打开  
     CFile fDbf;  
};

#endif //_MAP_RECORDSET_H_


/************************************************************
  
  MapRender
  
***********************************************************/

#ifndef _MAP_RENDER_H_
#define _MAP_RENDER_H_

#define MAX_CLASSNUM     //最大分类数

typedef enum 
{
	SIMPLE_RENDER,
    UNIQUE_RENDER
} RENDERTYPE;

typedef struct simpleRender
{
	COLORREF FillColor;
    COLORREF OutlineColor;
	int      iIndex;
}SIMPLERENDER;

typedef struct renderInfo
{
	CString   csValue;
	COLORREF  clr;

}RENDERINFO;

class CMapRender 
{
public:
   CMapRender();
   ~CMapRender();
public:
   void  Add(RENDERINFO& rInfo );
   void  RemoveByValue(CString& csValue);
   void  RemoveByIndex(int iIndex );
   int   GetCount();
   RENDERINFO* GetByValue(CString& csValue);
   RENDERINFO* GetByIndex(int iIndex);
   void Clear();
   void SetRenderType(int m_Type) { m_RenderType = m_Type; };
   int  GetRenderType() { return m_RenderType; };
   void SetSimpleRender( SIMPLERENDER& simpleRender );
   void GetSimpleRender( SIMPLERENDER& simpleRender );
   void SetFieldIndex(int iIndex);
   int  GetFieldIndex();
   void Clone(CMapRender *pRender);
private:  
   SIMPLERENDER  m_SimpleRender; 
   int m_RenderType;
   int m_FieldIndex;
private:
	CArray<RENDERINFO*,RENDERINFO*> m_Render;

};
#endif //_MAP_RENDER_H_


/************************************************************
  
  shpFile
  
***********************************************************/


#ifndef _SHPFILE_H_
#define _SHPFILE_H_

#define FILE_VERSION 1000
#define MAX_BUFFER_SIZE 32768 /*32K*/

/*系统错误常量*/

#define FILE_READERR    -1 
#define FILE_CODEERR    -2 
#define FILE_VERSIONERR -3

//shp 类型
#define NULLSHP       0
#define POINT         1
#define POLYLINE      3
#define POLYGON       5
#define MULTIPOINT    8
#define POINTZ        11
#define POLYLINEZ     13
#define POLYGONZ      15
#define MULTIPOINTZ   18
#define POINTM        21
#define POLYLINEM     23
#define POLYGONM      25
#define MULTIPOINTM   28
#define MULTIPATCH    31 

#pragma pack(4) 

typedef struct  {
	int iFileCode;      //文件标识
	int iReserved[5];   //保留字节
	int iFileLength;    //文件长度
	int iVersion;       //版本号
	int iShpType;       //文件类型
	double dbXMin;
	double dbYMin;
	double dbXMax;
	double dbYMax;
	double dbZMin;
	double dbZMax;
	double dbMMin;
	double dbMMax;
} SHPHEADER;

typedef struct shpRecordHeader {
	int iRecordNum;      //记录数
	int iContentLength;  //记录内容长度
} SHPRECORDHEADER;

typedef struct shpPoint
{
	double dbX;
	double dbY;
} SHPPOINT;

typedef struct shpPtContent {
	int iShpType;
	double dbX;
	double dbY;
} SHPPTCONTENT;


typedef struct shxRecord {
	int iOffset;
	int iContentLength;
} SHXRECORD;

typedef struct shpInfo {
	int ishpType;          //shp 类型
	shpPoint Box[2];       //最大矩形范围 
	int iNumParts;         //分段数
	int iNumPoints;        //shp总共顶点数  
} SHPINFO;

class CMapLayer;

class CShpFile 
{
public:
    CShpFile();
	virtual ~CShpFile();

public:
	BOOL ReadShp(CString& csFileName);
	CMapRectangle GetExtent();
	void SetExtent(CMapRectangle& extent );
	int GetShpType( ) {return m_shpType;};
	void SetShpType(int& iShpType);
	int SearchShape(CMapPoint& pt );
protected:
	BOOL ReadRecord();
	BOOL GetRecordHeader(SHPRECORDHEADER& varHeader );
	BOOL GetShpInfo(SHPINFO& varInfo);
	int  ReadPoint(char* pszBuffer,int& iLength,BOOL& bEof); 
	int SetRecordPos( int iRecord );
	BOOL ReadDBF(CString& csFileName);
	osg::ref_ptr<osg::Group> DrawShp(int mode);
private:
	void GetCurrLinePointElevation(Geodesy *geod, int mode);
	void GetCurrPolygonPointElevation(Geodesy *geod, int mode);
	void DrawPoint(int mode);
	void DrawPointElement(CMapPoint *pPoint, SimPoint *sp);
	osg::ref_ptr<osg::Group> DrawPolygon(int mode);
	osg::ref_ptr<osg::Group> DrawPolygonElement(int mode, CMapPolygon* pPolygon);
	osg::ref_ptr<osg::Group> DrawPLine(int mode);
	osg::ref_ptr<osg::Group> DrawPLineElement(int mode, CMapLine* pPline, float roadWidth, int type);
private:
	BOOL bShpOpen,bShxOpen;
	BOOL m_bBigEndian; 
	CFile fShp,fShx;
	CMapRectangle m_Extent;
	CMapRecordSet *pRecordSet;
private:
	int m_shpFileLength;
	int m_shxFileLength;
	int m_iRecordCount;
	int m_shpType;
	CMapRectangle m_CurMapExtent;
	CObList m_ObList;
	bool isLake;
	bool islandInOcean;
	bool isAirport;
	double m_elevOfLake;
	double m_elevOfAirport;
private:
	friend class CMapLayer;
};
#endif //_SHPFILE_H_

/************************************************************
  
  MapLayer
  
***********************************************************/

#ifndef _MAP_LAYER_H_
#define _MAP_LAYER_H_

class CMapLayer {

public:
	CMapLayer();
	~CMapLayer();

//Attribute
public:
    CMapRectangle GetExtent();
	BOOL ReadShp(CString& csFileName);
	long GetShpType();
	BOOL GetVisible();
	BOOL LoadData(CString& csFileName);
	void SetVisible(bool bVisible);
	void SetLayerName(CString& csLayerName );
	CString GetLayerName();
    int SearchShape(CMapPoint& pt );
	void SetRender(CMapRender* pRender ); 
	osg::ref_ptr<osg::Group> DrawLayer(int mode);
	CMapRecordSet* GetRecordSet() { return m_shpFile.pRecordSet;};
    CMapRender* GetRender() { return m_pRender; } ;
protected:
	void SetExtent(CMapRectangle& extent );
	void SetShpType(int lShpType );

public:
	
	//CMapRecordSet *pRecordSet;
private:  
	BOOL m_bVisible;
	BOOL m_Valid;
	CMapRectangle m_Extent;
	CShpFile m_shpFile;
	CString m_csLayerName;
	CMapRender* m_pRender;
};	

#endif //_MAP_LAYER_H_ 

/************************************************************
 	
 	MapLayers
 	
***********************************************************/

#ifndef _MAP_LAYERS_H_
#define _MAP_LAYERS_H_

class CMapLayers {
public:
	CMapLayers();
	~CMapLayers();

// Attributes
public:
	short GetCount();
	void SetCount(short);
	void GetAllExtent(CMapRectangle& rc);

// Operations
public:
	CMapLayer* GetAt(short sIndex);
	BOOL Add(CMapLayer* pMapLayer);
	void Remove(short sIndex);
	void Clear();
	void MoveTo(short fromIndex, short toIndex);
	void MoveToTop(short sIndex);
	void MoveToBottom(short sIndex);
private:
	CArray<CMapLayer*,CMapLayer*> m_Layers; //存储图层集合对象
};

#endif //_MAP_LAYERS_H_

/************************************************************

 	Interface

***********************************************************/

int LoadShpAndDrawIt(char *shpfn);
int LoadShpAndCalculateMM(char *shpfn);
