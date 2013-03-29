#pragma once


// CNumEdit

class CNumEdit : public CEdit
{
	DECLARE_DYNAMIC(CNumEdit)

public:
	CNumEdit();
	virtual ~CNumEdit();
	
	BOOL CheckNumber(UINT nChar,UINT nRepCnt,UINT nFlags);
	BOOL CheckOneMinus(UINT nChar,UINT nRepCnt,UINT nFlags);
	BOOL CheckOneDot(UINT nChar,UINT nRepCnt,UINT nFlags);
	int  GetCaretXPos();
    double   GetDouble();
	long GetLong();
	void SetWindowText( CString str );
    CString m_editstr;
	int m_iAfterDotLen;
protected:
	DECLARE_MESSAGE_MAP()
	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnKillFocus(CWnd* pNewWnd);
};


