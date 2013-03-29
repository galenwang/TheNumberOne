#include "dijingzhizuo.h"
#if   !defined(AFX_MYTABCTRL_H__E5566EC0_70E7_4D7C_A413_967F5DEB779A__INCLUDED_) 
#define   AFX_MYTABCTRL_H__E5566EC0_70E7_4D7C_A413_967F5DEB779A__INCLUDED_ 

#if   _MSC_VER   >   1000 
#pragma   once 
#endif   //   _MSC_VER   >   1000 
//   MyTabCtrl.h   :   header   file 
// 

///////////////////////////////////////////////////////////////////////////// 
//   CMyTabCtrl   window 

class   CMyTabCtrl   :   public   CTabCtrl 
{ 
	//   Construction 
public: 
	CMyTabCtrl(); 

	//   Attributes 
public: 

	//   Operations 
public: 

	//   Overrides 
	//   ClassWizard   generated   virtual   function   overrides 
	//{{AFX_VIRTUAL(CMyTabCtrl) 
public: 
	//}}AFX_VIRTUAL 

	//   Implementation 
public: 
	virtual   ~CMyTabCtrl(); 
	CImageList* pImageList;

	//   Generated   message   map   functions 
protected: 
	//{{AFX_MSG(CMyTabCtrl) 
	afx_msg   void   OnPaint(); 
	afx_msg   BOOL   OnEraseBkgnd(CDC*   pDC); 
	//}}AFX_MSG 

	DECLARE_MESSAGE_MAP() 
	CDC   m_memDC; 
	CBitmap   m_memBp; 
}; 

///////////////////////////////////////////////////////////////////////////// 

//{{AFX_INSERT_LOCATION}} 
//   Microsoft   Visual   C++   will   insert   additional   declarations   immediately   before   the   previous   line. 

#endif   //   !defined(AFX_MYTABCTRL_H__E5566EC0_70E7_4D7C_A413_967F5DEB779A__INCLUDED_) 
