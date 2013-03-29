// NumEdit.cpp : ʵ���ļ�
//

#include "stdafx.h"
//#include "instructor operator station.h"
#include "NumEdit.h"


// CNumEdit

IMPLEMENT_DYNAMIC(CNumEdit, CEdit)

CNumEdit::CNumEdit()
{
 m_iAfterDotLen = 6;
 m_editstr = _T("0");
}

CNumEdit::~CNumEdit()
{
}

BEGIN_MESSAGE_MAP(CNumEdit, CEdit)
 //{{AFX_MSG_MAP(CNumberEdit)
 ON_WM_CHAR()
 ON_WM_KILLFOCUS()

// ON_CONTROL_REFLECT(EN_KILLFOCUS,OnKillFocus)
 //}}AFX_MSG_MAP
END_MESSAGE_MAP()

void CNumEdit::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
 // �޸���Ϣ��Ӧ
 if(nChar == 8)
 {
  //�˸�
  CEdit::OnChar(nChar, nRepCnt, nFlags);
  return;
 }
 BOOL bChange = FALSE;
 GetWindowText(m_editstr);
 if(CheckNumber(nChar,nRepCnt,nFlags))
 {
  //�������Ϊ����
  bChange = TRUE;
 }
 else if(CheckOneMinus(nChar,nRepCnt,nFlags))
 {
  //���ֻ��һ������,��ֻ���ǵ�һ���ַ�
  bChange = TRUE;
 }
 else if(CheckOneDot(nChar,nRepCnt,nFlags))
 {
  //���ֻ��һ��С����
  bChange = TRUE;
 }
 if(bChange)
 {
  CEdit::OnChar(nChar, nRepCnt, nFlags);
 }
}

void CNumEdit::OnKillFocus(CWnd* pNewWnd) 
{
 CEdit::OnKillFocus(pNewWnd);

 // �޸���Ϣ��Ӧ
 //ʧȥ����ʱ��С�������õ�С��������ݳ��Ȳ���
//*
GetWindowText(m_editstr);
 if(m_editstr.IsEmpty())  
 {
  //δ������������Ϊ0
  m_editstr = "0";
 }
 else if(m_editstr.GetLength() == 1 && m_editstr[0] == TCHAR('-'))
 {
  //ֻ������һ������
  m_editstr = "0";
 }
 else
 {
  int iDotPos = m_editstr.Find(TCHAR('.'));
  if(iDotPos <0)
  {
   //û���ҵ�С����
   SetWindowText(m_editstr);
   return;
  }
  if(iDotPos >0)
  {
   //С���������λ��
   int iLen = m_editstr.GetLength() - 1 - iDotPos;
   if(iLen >= m_iAfterDotLen)
   {
    //����λ�������趨
    return;
   }
   if (iLen == 0)
   {
    //С�����û����������0����
    m_editstr+="0";
   }
  }
 }
 SetWindowText(m_editstr);
// */
}

BOOL CNumEdit::CheckNumber(UINT nChar,UINT nRepCnt,UINT nFlags)
{
 if(::isdigit(nChar)==0)
 {
  //�Ƿ�����
  return FALSE;
 }
 //С����λ��
 int iDotPos = m_editstr.Find(TCHAR('.'));
 if(iDotPos >= 0)
 {
  //С���ݵ�����ݳ���
  int iLen = m_editstr.GetLength() - 1 - iDotPos;
  if( (GetCaretXPos() >= iDotPos) && (iLen >= m_iAfterDotLen))
  {
   //�������õĳ���
   return FALSE;
  }
 }
 return TRUE;
}


BOOL CNumEdit::CheckOneMinus(UINT nChar,UINT nRepCnt,UINT nFlags)
{
 if(nChar != '-')
 {
  //����'-'
  return FALSE;
 }
 if(GetCaretXPos() != 0)
 {
  //���ڵ�һ��λ��
  return FALSE;
 }
 if(!m_editstr.IsEmpty() && m_editstr.GetAt(0) == TCHAR('-'))
 {
  return FALSE;
 }
 return TRUE;
}

BOOL CNumEdit::CheckOneDot(UINT nChar,UINT nRepCnt,UINT nFlags)
{
 if (m_iAfterDotLen == 0)
 {
  //С����󲻼�����(����Ϊ����)
  return FALSE;
 }
 if(nChar != '.')
 {
  //����С����
  return FALSE;
 }
 if(m_editstr.Find(TCHAR('.')) >=0)
 {
  //����С����
  return FALSE;
 }
 int iPos = GetCaretXPos();
 if(iPos == 0)
 {
  //��һ���ַ�����С���㣡
  return FALSE;
 }
 else if(iPos==1 && m_editstr[0] == TCHAR('-'))
 {
  //��һ���ַ���'-',�ڶ�����С����
  return FALSE;
 }
 return TRUE;
}

int CNumEdit::GetCaretXPos()
{
 CPoint p = GetCaretPos();
 return (p.x - p.y)/6;
}

double CNumEdit::GetDouble()
{
 GetWindowText(m_editstr);
 return atof((const char *)(LPCTSTR)m_editstr.GetBuffer());
}

long CNumEdit::GetLong()
{
 GetWindowText(m_editstr);
 return atol((const char *)(LPCTSTR)m_editstr.GetBuffer());
}

void CNumEdit::SetWindowText( CString str )
{
 m_editstr = str;
 char nChar;
 UINT nRepCnt = 0;
 UINT nFlags = 0;
 BOOL bChange = FALSE;
 for (int nstr = 0;nstr < str.GetLength() - 1;++nstr)
 {
  //ȡ���ַ�
  nChar = str.GetAt(nstr);
  if(!CheckNumber(nChar,nRepCnt,nFlags) && (!CheckOneMinus

(nChar,nRepCnt,nFlags))
   && (!CheckOneDot(nChar,nRepCnt,nFlags)))
  {
   //����ʼ���ı�����
   bChange = TRUE;
  }
 }
 if(!bChange)
 {
  int iDotPos = m_editstr.Find(TCHAR('.'));
  if(iDotPos >0)
  {
   //С���������λ��
   int iLen = m_editstr.GetLength() - 1 - iDotPos;
   if (iLen == 0)
   {
    //С�����û����������0����
    m_editstr+="0";
   }
  }
  CEdit::SetWindowText(m_editstr);
 }
}


// CNumEdit ��Ϣ�������


