// NumEdit.cpp : 实现文件
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
 // 修改消息响应
 if(nChar == 8)
 {
  //退格
  CEdit::OnChar(nChar, nRepCnt, nFlags);
  return;
 }
 BOOL bChange = FALSE;
 GetWindowText(m_editstr);
 if(CheckNumber(nChar,nRepCnt,nFlags))
 {
  //检查输入为数字
  bChange = TRUE;
 }
 else if(CheckOneMinus(nChar,nRepCnt,nFlags))
 {
  //检查只有一个负号,且只能是第一个字符
  bChange = TRUE;
 }
 else if(CheckOneDot(nChar,nRepCnt,nFlags))
 {
  //检查只有一个小数点
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

 // 修改消息响应
 //失去焦点时，小数按设置的小数点后数据长度补齐
//*
GetWindowText(m_editstr);
 if(m_editstr.IsEmpty())  
 {
  //未输入数据则设为0
  m_editstr = "0";
 }
 else if(m_editstr.GetLength() == 1 && m_editstr[0] == TCHAR('-'))
 {
  //只输入了一个负号
  m_editstr = "0";
 }
 else
 {
  int iDotPos = m_editstr.Find(TCHAR('.'));
  if(iDotPos <0)
  {
   //没有找到小数点
   SetWindowText(m_editstr);
   return;
  }
  if(iDotPos >0)
  {
   //小数点后已有位数
   int iLen = m_editstr.GetLength() - 1 - iDotPos;
   if(iLen >= m_iAfterDotLen)
   {
    //已有位数超过设定
    return;
   }
   if (iLen == 0)
   {
    //小数点后没有数据则用0补齐
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
  //是否数字
  return FALSE;
 }
 //小数点位置
 int iDotPos = m_editstr.Find(TCHAR('.'));
 if(iDotPos >= 0)
 {
  //小数据点后数据长度
  int iLen = m_editstr.GetLength() - 1 - iDotPos;
  if( (GetCaretXPos() >= iDotPos) && (iLen >= m_iAfterDotLen))
  {
   //超过设置的长度
   return FALSE;
  }
 }
 return TRUE;
}


BOOL CNumEdit::CheckOneMinus(UINT nChar,UINT nRepCnt,UINT nFlags)
{
 if(nChar != '-')
 {
  //不是'-'
  return FALSE;
 }
 if(GetCaretXPos() != 0)
 {
  //不在第一个位置
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
  //小数点后不加数据(限制为整数)
  return FALSE;
 }
 if(nChar != '.')
 {
  //不是小数点
  return FALSE;
 }
 if(m_editstr.Find(TCHAR('.')) >=0)
 {
  //已有小数点
  return FALSE;
 }
 int iPos = GetCaretXPos();
 if(iPos == 0)
 {
  //第一个字符就是小数点！
  return FALSE;
 }
 else if(iPos==1 && m_editstr[0] == TCHAR('-'))
 {
  //第一个字符是'-',第二个是小数点
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
  //取得字符
  nChar = str.GetAt(nstr);
  if(!CheckNumber(nChar,nRepCnt,nFlags) && (!CheckOneMinus

(nChar,nRepCnt,nFlags))
   && (!CheckOneDot(nChar,nRepCnt,nFlags)))
  {
   //检查初始化文本数据
   bChange = TRUE;
  }
 }
 if(!bChange)
 {
  int iDotPos = m_editstr.Find(TCHAR('.'));
  if(iDotPos >0)
  {
   //小数点后已有位数
   int iLen = m_editstr.GetLength() - 1 - iDotPos;
   if (iLen == 0)
   {
    //小数点后没有数据则用0补齐
    m_editstr+="0";
   }
  }
  CEdit::SetWindowText(m_editstr);
 }
}


// CNumEdit 消息处理程序


