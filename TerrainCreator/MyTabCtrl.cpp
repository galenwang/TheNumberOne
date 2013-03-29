//   MyTabCtrl.cpp   :   implementation   file 
// 

#include "stdafx.h" 

#include   "MyTabCtrl.h" 

#ifdef   _DEBUG 
#define   new   DEBUG_NEW 
#undef   THIS_FILE 
static   char   THIS_FILE[]   =   __FILE__; 
#endif 

///////////////////////////////////////////////////////////////////////////// 
//   CMyTabCtrl 

CMyTabCtrl::CMyTabCtrl() 
{ 
} 

CMyTabCtrl::~CMyTabCtrl() 
{
	if(m_memDC.GetSafeHdc()) 
	m_memDC.DeleteDC(); 
	if(m_memBp.GetSafeHandle()) 
	m_memBp.DeleteObject(); 
}


BEGIN_MESSAGE_MAP(CMyTabCtrl,   CTabCtrl) 
	//{{AFX_MSG_MAP(CMyTabCtrl) 
	ON_WM_PAINT() 
	ON_WM_ERASEBKGND() 
	//}}AFX_MSG_MAP 
END_MESSAGE_MAP() 

///////////////////////////////////////////////////////////////////////////// 
//   CMyTabCtrl   message   handlers 

void   CMyTabCtrl::OnPaint()   
{ 
	// CTabCtrl::Default(); 
	static   BOOL   bFirst   =   TRUE; 
	static   CRect   rcTabInParent; 

	CPaintDC   dc(this);   //   device   context   for   painting 
	CRect   rc,   rcTab; 
	GetItemRect(0,   &rc); 
	GetWindowRect(&rcTab); 
	ScreenToClient(&rcTab); 
	rc.right   =   rcTab.right; 
	rc.left   =   rc.left   -   2; 
//	rc.top = rc.top -100;

	if(bFirst)//这里把父窗口所画的保存起来 
	{ 
		bFirst   =   FALSE; 
		rcTabInParent   =   rc; 
		rcTabInParent.left   =   rcTabInParent.top   =   0; 
		CWnd*   pParent   =   GetParent(); 
		CDC*   pParentDC   =   pParent-> GetDC(); 
		m_memBp.CreateCompatibleBitmap(&dc,   rcTabInParent.Width(),   rcTabInParent.Height()); 
		m_memDC.CreateCompatibleDC(&dc); 
		m_memDC.SelectObject(m_memBp.GetSafeHandle()); 

		ClientToScreen(&rcTabInParent); 
		pParent-> ScreenToClient(&rcTabInParent); 

	//	m_memDC.StretchBlt(0,   0,   rcTabInParent.Width(),   rcTabInParent.Height(),   pParentDC,   rcTabInParent.left,   rcTabInParent.top,   rcTabInParent.Width(),   rcTabInParent.Height(),SRCCOPY); 
		pParent-> ReleaseDC(pParentDC); 
	} 

	//dc.BitBlt(0,  0,   rcTabInParent.Width(),   rcTabInParent.Height(),   &m_memDC,   0,   0,   SRCCOPY);//把所保存的父窗口所画的画在此控件上 

	int   iSel   =   GetCurSel(); 
 	for(int   i   =   0;   i   <   GetItemCount();   i++)//这里画控件的Item及相应的状态 
 	{ 
 		GetItemRect(i,   &rc); 
 		TCITEM   Item; 
 		Item.mask   =   TCIF_TEXT; 
		Item.pszText   =   new   TCHAR[100]; 
 		Item.cchTextMax   =   100; 
 		ZeroMemory(Item.pszText,   sizeof(TCHAR)   *   100); 
 		GetItem(i,   &Item);
 	/*	if(iSel   ==   i) 
 		{ 
 			dc.DrawEdge(&rc,   EDGE_BUMP,   BF_TOPLEFT   |   BF_TOPRIGHT); 
 		} 
 		else 
 		{ 
 			//rc.top   -=   3;
			rc.left-=0.5;
 			dc.DrawEdge(&rc,   EDGE_BUMP,   BF_RECT); 
 		} */
//  		int   iOldMode   =   dc.SetBkMode(TRANSPARENT); 
//  		HGDIOBJ   hOldObj   =   dc.SelectObject(GetParent()-> GetFont()); 
//  		dc.DrawText(Item.pszText,   &rc,   DT_CENTER   |   DT_VCENTER   |   DT_SINGLELINE   |   DT_END_ELLIPSIS); 
//  		dc.SetBkMode(iOldMode); 
//  		dc.SelectObject(hOldObj); 
 		delete[]   Item.pszText; 
	} 

	//   TODO:   Add   your   message   handler   code   here 
 	

	//   Do   not   call   CTabCtrl::OnPaint()   for   painting   messages */
} 

BOOL   CMyTabCtrl::OnEraseBkgnd(CDC*   pDC)   
{ 
	//   TODO:   Add   your   message   handler   code   here   and/or   call   default 

	return   TRUE; 
	// return   CTabCtrl::OnEraseBkgnd(pDC); 
} 
