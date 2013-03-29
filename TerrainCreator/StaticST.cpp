#include "stdafx.h"
#include "StaticST.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CStaticST::CStaticST()
{
	FreeResources(FALSE);

	// Default colors
	m_crTextColor = ::GetSysColor(COLOR_BTNTEXT);
	m_crBkColor = ::GetSysColor(COLOR_3DFACE);

	// No messages to sent on clicks
	::ZeroMemory(&m_csClicks, sizeof(m_csClicks));

	// Do not draw as a transparent button
	m_bDrawTransparent = FALSE;
	m_pbmpOldBk = NULL;
}

CStaticST::~CStaticST()
{
	// Restore old bitmap (if any)
	if (m_dcBk.m_hDC && m_pbmpOldBk)
	{
		m_dcBk.SelectObject(m_pbmpOldBk);
	} // if

	FreeResources();
}

BEGIN_MESSAGE_MAP(CStaticST, CStatic)
	//{{AFX_MSG_MAP(CStaticST)
	ON_WM_PAINT()
	ON_WM_NCPAINT()
	ON_WM_LBUTTONUP()
	ON_WM_LBUTTONDBLCLK()
	ON_WM_CTLCOLOR_REFLECT()
	ON_WM_ERASEBKGND()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

void CStaticST::FreeResources(BOOL bCheckForNULL)
{
	if (bCheckForNULL)
	{
		// Destroy the icon (if any)
		// Note: the following line MUST be here! even if
		// BoundChecker says it's unnecessary!
		if (m_hIcon != NULL) ::DestroyIcon(m_hIcon);
	} // if

	// No icon
	m_hIcon = NULL;
	m_cxIcon = 0;
	m_cyIcon = 0;
} // End of FreeResources

void CStaticST::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
	
	DrawItem(&dc);
} // End of OnPaint

void CStaticST::OnNcPaint() 
{
	// Do not call CStatic::OnNcPaint() for painting messages
} // End of OnNcPaint

BOOL CStaticST::OnEraseBkgnd(CDC* pDC) 
{
	return FALSE;
} // End of OnEraseBkgnd

HBRUSH CStaticST::CtlColor(CDC* pDC, UINT nCtlColor) 
{
	return (HBRUSH)::GetStockObject(NULL_BRUSH); 
} // End of CtlColor

void CStaticST::OnLButtonUp(UINT nFlags, CPoint point) 
{
	CStatic::OnLButtonUp(nFlags, point);

	// Send message to window
	if (m_csClicks[STATICST_CLICK].hWnd != NULL)
	{
		::PostMessage(m_csClicks[STATICST_CLICK].hWnd, m_csClicks[STATICST_CLICK].nMsg, 0, 0);
	} // if
} // End of OnLButtonUp

void CStaticST::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	CStatic::OnLButtonDblClk(nFlags, point);

	// Send message to window
	if (m_csClicks[STATICST_DBLCLICK].hWnd != NULL)
	{
		::PostMessage(m_csClicks[STATICST_DBLCLICK].hWnd, m_csClicks[STATICST_DBLCLICK].nMsg, 0, 0);
	} // if
} // End of OnLButtonDblClk

void CStaticST::DrawItem(CDC* pDC)
{
	CString sText;
	CRect rCtrl;

	GetClientRect(rCtrl);

	DWORD dwStyle	= GetStyle();
	DWORD dwExStyle	= GetExStyle();

	// Prepare draw... paint control background

	// Draw transparent?
	if (m_bDrawTransparent)
		PaintBk(pDC);
	else
		OnDrawBackground(pDC, &rCtrl);

	//if (dwStyle & SS_SUNKEN)			rCtrl.InflateRect(1, 1);
	//if (dwExStyle & WS_EX_CLIENTEDGE)	rCtrl.InflateRect(2, 2);

	// SS_RIGHT			Right aligned
	// SS_CENTER		Center aligned
	// SS_CENTERIMAGE	Center aligned (vertically)
	// SS_SUNKEN		Sunken
	// WS_EX_CLIENTEDGE	Border
	if (dwStyle & SS_SUNKEN)
	{
		pDC->Draw3dRect(rCtrl, ::GetSysColor(COLOR_BTNSHADOW), ::GetSysColor(COLOR_BTNHILIGHT));
		rCtrl.DeflateRect(1, 1);
	} // if

	if (dwExStyle & WS_EX_CLIENTEDGE)
	{
		pDC->Draw3dRect(rCtrl, RGB(0,0,0), RGB(0,0,0));
		rCtrl.DeflateRect(1, 1);
	} // if

	GetWindowText(sText);

	// Draw icon
	if (m_hIcon)	DrawTheIcon(pDC, &rCtrl, (BOOL)(sText.IsEmpty() == 0));

	// Draw text
	DrawTheText(pDC, &rCtrl, sText, dwStyle, dwExStyle);
} // End of DrawItem

void CStaticST::DrawTheIcon(CDC* pDC, CRect* rpCtrl, BOOL bText)
{
	int		nX = 0, nY = 0;
	CRect	rIcon;

	rIcon.CopyRect(rpCtrl);

	// If there is text
	if (bText == TRUE)
	{
		rIcon.right = rIcon.left + (4+m_cxIcon+4);
		rpCtrl->left = rIcon.right;
	}

	nX = rIcon.left;
	nY = rIcon.top;

	nX += (rIcon.Width() - m_cxIcon) / 2;
	//if (nX < 0) nX = 0;
	nY += (rIcon.Height() - m_cyIcon) / 2;
	//if (nY < 0) nY = 0;

	pDC->DrawState(	CPoint(nX, nY),
					CSize(m_cxIcon, m_cyIcon),
					m_hIcon,
					DSS_NORMAL,
					(CBrush*)NULL);
} // End of DrawTheIcon

void CStaticST::DrawTheText(CDC* pDC, CRect* rpCtrl, CString& sText, DWORD dwStyle, DWORD dwExStyle)
{
	CRect	rText;
	int		nX = 0, nY = 0;
	UINT	nFormat = DT_EXPANDTABS;

	// Transparent background
	pDC->SetBkMode(TRANSPARENT);

	pDC->SetTextColor(m_crTextColor);

	// Set same font of the window
	CFont* oldFont = pDC->SelectObject(GetFont());

	if ((dwStyle & SS_NOPREFIX) == SS_NOPREFIX)	nFormat |= DT_NOPREFIX;
	if ((dwStyle & SS_LEFTNOWORDWRAP) != SS_LEFTNOWORDWRAP) nFormat |= DT_WORDBREAK;

	rText = rpCtrl;
	pDC->DrawText(sText, -1, rText, nFormat | DT_CALCRECT);

	// Center horizontally
	if (dwStyle & SS_CENTER)
	{
		nX += (rpCtrl->Width() - rText.Width())/2;
	} // if
	
	// Center vertically
	if (dwStyle & SS_CENTERIMAGE)
	{
		nY += (rpCtrl->Height() - rText.Height())/2;
	} // if

	rText.OffsetRect(nX, nY);

	pDC->DrawText(sText, -1, rText, nFormat);

	// ONLY FOR DEBUG 
	//CBrush brBtnShadow(RGB(255, 0, 0));
	//pDC->FrameRect(&rText, &brBtnShadow);

	// Restore old font
	pDC->SelectObject(oldFont);
} // End of DrawTheText

// This function sets foreground color of the control.
//
// Parameters:
//		[IN]	crTextColor
//				A COLORREF value indicating the foreground color
//
//		[IN]	bRepaint
//				If TRUE the control will be repainted.
//
void CStaticST::SetTextColor(COLORREF crTextColor, BOOL bRepaint)
{
	m_crTextColor = crTextColor;

	// Repaint control
	if (bRepaint)	Invalidate();
} // End of SetTextColor

// This function sets background color of the control.
//
// Parameters:
//		[IN]	crBkColor
//				A COLORREF value indicating the background color
//
//		[IN]	bRepaint
//				If TRUE the control will be repainted.
//
void CStaticST::SetBkColor(COLORREF crBkColor, BOOL bRepaint)
{
	m_crBkColor = crBkColor;

	// Repaint control
	if (bRepaint)	Invalidate();
} // End of SetBkColor

// This function sets foreground and background colors of the control.
//
// Parameters:
//		[IN]	crTextColor
//				A COLORREF value indicating the foreground color
//
//		[IN]	crBkColor
//				A COLORREF value indicating the background color
//
//		[IN]	bRepaint
//				If TRUE the control will be repainted.
//
void CStaticST::SetColors(COLORREF crTextColor, COLORREF crBkColor, BOOL bRepaint)
{
	SetTextColor(crTextColor, FALSE);
	SetBkColor(crBkColor, bRepaint);
} // End of SetColors

// This function sets the icon to be displayed.
// Any previous icon will be removed.
//
// Parameters:
//		[IN]	nIcon
//				A Windows icon resource ID
//		[IN]	bRepaint
//				If TRUE the control will be immediately repainted
//		[IN]	hInstance
//				Handle of the instance that contains the icon.
//				If NULL the icon will be loaded from the .EXE resources
//
// Return value:
//		STATICST_OK
//			Function executed successfully.
//		STATICST_INVALIDRESOURCE
//			Failed loading the specified resource.
//
DWORD CStaticST::SetIcon(int nIcon, BOOL bRepaint, HINSTANCE hInstance)
{
	HINSTANCE	hInstResource;
	HICON		hIcon;

	if (hInstance == NULL)
		hInstResource = AfxFindResourceHandle(MAKEINTRESOURCE(nIcon), RT_GROUP_ICON);
	else
		hInstResource = hInstance;

	// Load icon from resource
	hIcon = (HICON)::LoadImage(hInstResource, MAKEINTRESOURCE(nIcon), IMAGE_ICON, 0, 0, 0);

	return SetIcon(hIcon, bRepaint);
} // End of SetIcon

// This function sets the icon to be displayed.
// Any previous icon will be removed.
//
// Parameters:
//		[IN]	hIcon
//				Handle fo the icon to show.
//				Pass NULL to remove any icon from the control.
//		[IN]	bRepaint
//				If TRUE the control will be immediately repainted
//
// Return value:
//		STATICST_OK
//			Function executed successfully.
//		STATICST_INVALIDRESOURCE
//			Failed loading the specified resource.
//
DWORD CStaticST::SetIcon(HICON hIcon, BOOL bRepaint)
{
	BOOL		bRetValue;
	ICONINFO	ii;

	// Free any loaded resource
	FreeResources();

	if (hIcon)
	{
		m_hIcon = hIcon;

		// Get icon dimension
		ZeroMemory(&ii, sizeof(ICONINFO));
		bRetValue = ::GetIconInfo(m_hIcon, &ii);
		if (bRetValue == FALSE)
		{
			FreeResources();
			return STATICST_INVALIDRESOURCE;
		} // if

		m_cxIcon = (BYTE)(ii.xHotspot * 2);
		m_cyIcon = (BYTE)(ii.yHotspot * 2);
		::DeleteObject(ii.hbmMask);
		::DeleteObject(ii.hbmColor);
	} // if

	// Repaint control
	if (bRepaint)	Invalidate();

	return STATICST_OK;
} // End of SetIcon

// Sets bold effect for control's text
//
// Parameters:
//		[IN]	nBold
//				If TRUE the control's text will have a bold effect
//
//		[IN]	bRepaint
//				If TRUE the control will be immediately repainted
//
// Return value:	none
//
void CStaticST::SetBold(BOOL bBold, BOOL bRepaint)
{
	LOGFONT logfont;
	CFont*	pFont = NULL;

	pFont = GetFont();
	if (pFont)
	{
		if (pFont->GetLogFont(&logfont) != 0)
		{
			logfont.lfWeight = (bBold == TRUE) ? FW_BOLD:FW_NORMAL;
			m_Font.DeleteObject();
			m_Font.CreateFontIndirect(&logfont);
			SetFont(&m_Font, bRepaint);
		} // if
	} // if
} // End of SetBold

void CStaticST::InitToolTip()
{
	if (m_ToolTip.m_hWnd == NULL)
	{
		// Create ToolTip control
		m_ToolTip.Create(this);
		// Create inactive
		m_ToolTip.Activate(FALSE);
		// Enable multiline
		m_ToolTip.SendMessage(TTM_SETMAXTIPWIDTH, 0, 400);
	} // if
} // End of InitToolTip

void CStaticST::ActivateTooltip(BOOL bActivate)
{
	// If there is no tooltip then do nothing
	if (m_ToolTip.GetToolCount() == 0) return;

	// Activate tooltip
	m_ToolTip.Activate(bActivate);
} // End of EnableTooltip

void CStaticST::SetTooltipText(LPCTSTR lpszText, BOOL bActivate)
{
	// We cannot accept NULL pointer
	if (lpszText == NULL) return;

	// Initialize ToolTip
	InitToolTip();

	// If there is no tooltip defined then add it
	if (m_ToolTip.GetToolCount() == 0)
	{
		CRect rectBtn; 
		GetClientRect(rectBtn);
		m_ToolTip.AddTool(this, lpszText, rectBtn, 1);
	} // if

	// Set text for tooltip
	m_ToolTip.UpdateTipText(lpszText, this, 1);
	m_ToolTip.Activate(bActivate);
} // End of SetTooltipText

void CStaticST::SetTooltipText(int nId, BOOL bActivate)
{
	CString sText;

	// load string resource
	sText.LoadString(nId);
	// If string resource is not empty
	if (sText.IsEmpty() == FALSE) SetTooltipText((LPCTSTR)sText, bActivate);
} // End of SetTooltipText

BOOL CStaticST::PreTranslateMessage(MSG* pMsg) 
{
	InitToolTip();
	m_ToolTip.RelayEvent(pMsg);
	
	return CStatic::PreTranslateMessage(pMsg);
} // End of PreTranslateMessage

// This function assigns a message that will be sent to a window when the user
// clicks on the control.
//
// Parameters:
//		[IN]	byClickType
//				Type of the click for wich assign the message. This is a
//				zero-based value and can be one of the following:
//				STATICST_CLICK		A single click on the control
//				STATICST_DBLCLICK	A double click on the control
//		[IN]	hWnd
//				Handle of the window that will be receive the message.
//				Set this value to NULL to disable the message for a
//				particular click.
//		[IN]	nMsg
//				Message to sent.
// Return value:
//		STATICST_OK
//			Function executed successfully.
//		STATICST_INVALIDCLICKTYPE
//			Invalid click type.
//
DWORD CStaticST::SetMessageClick(BYTE byClickType, HWND hWnd, UINT nMsg)
{
	// Check if valid click type
	if (byClickType >= STATICST_MAX_CLICKS)	return STATICST_INVALIDCLICKTYPE;

	m_csClicks[byClickType].hWnd = hWnd;
	m_csClicks[byClickType].nMsg = nMsg;

	return STATICST_OK;
} // End of SetMessageClick

void CStaticST::PaintBk(CDC* pDC)
{
	CClientDC clDC(GetParent());
	CRect rect;
	CRect rect1;

	GetClientRect(rect);

	GetWindowRect(rect1);
	GetParent()->ScreenToClient(rect1);

	if (m_dcBk.m_hDC == NULL)
	{
		m_dcBk.CreateCompatibleDC(&clDC);
		m_bmpBk.CreateCompatibleBitmap(&clDC, rect.Width(), rect.Height());
		m_pbmpOldBk = m_dcBk.SelectObject(&m_bmpBk);
		m_dcBk.BitBlt(0, 0, rect.Width(), rect.Height(), &clDC, rect1.left, rect1.top, SRCCOPY);
	} // if

	pDC->BitBlt(0, 0, rect.Width(), rect.Height(), &m_dcBk, 0, 0, SRCCOPY);
} // End of PaintBk

DWORD CStaticST::SetBk(CDC* pDC)
{
	if (m_bDrawTransparent && pDC)
	{
		// Restore old bitmap (if any)
		if (m_dcBk.m_hDC != NULL && m_pbmpOldBk != NULL)
		{
			m_dcBk.SelectObject(m_pbmpOldBk);
		} // if

		m_bmpBk.DeleteObject();
		m_dcBk.DeleteDC();

		CRect rect;
		CRect rect1;

		GetClientRect(rect);

		GetWindowRect(rect1);
		GetParent()->ScreenToClient(rect1);

		m_dcBk.CreateCompatibleDC(pDC);
		m_bmpBk.CreateCompatibleBitmap(pDC, rect.Width(), rect.Height());
		m_pbmpOldBk = m_dcBk.SelectObject(&m_bmpBk);
		m_dcBk.BitBlt(0, 0, rect.Width(), rect.Height(), pDC, rect1.left, rect1.top, SRCCOPY);

		return STATICST_OK;
	} // if

	return STATICST_BADPARAM;
} // End of SetBk

DWORD CStaticST::DrawTransparent(BOOL bRepaint)
{
	m_bDrawTransparent = TRUE;

	// Restore old bitmap (if any)
	if (m_dcBk.m_hDC != NULL && m_pbmpOldBk != NULL)
	{
		m_dcBk.SelectObject(m_pbmpOldBk);
	} // if

	m_bmpBk.DeleteObject();
	m_dcBk.DeleteDC();

	// Repaint the button
	if (bRepaint) Invalidate();

	return STATICST_OK;
} // End of DrawTransparent

DWORD CStaticST::OnDrawBackground(CDC* pDC, LPCRECT pRect)
{
	// Draw control background
	CBrush br(m_crBkColor);
	pDC->FillRect(pRect, &br);

	return STATICST_OK;
} // End of OnDrawBackground
