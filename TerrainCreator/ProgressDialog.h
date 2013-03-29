#pragma once


// CProgressDialog 对话框
class COperEnv;

class CProgressDialog : public CDialog
{
	DECLARE_DYNAMIC(CProgressDialog)

public:
	CProgressDialog(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CProgressDialog();
    void SetParentPointer(COperEnv* pParent);
// 对话框数据
	enum { IDD = IDD_PROGRESSDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnPaint();
	afx_msg void OnBnClickedOk();
	bool mark;
	int i,j;
	afx_msg void OnBnClickedButton1();
	int m_XSum,m_YSum,m_X,m_Y;
	virtual BOOL OnInitDialog();
protected:
	COperEnv*  theParent;
};
