// LatLon.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "LatLon.h"
#include "ProgressDialog.h"
#include "interface.h"
#include "travelmap.h"

extern float  mWoodCoverage;
extern float  mRoadWidthInShp;
extern float  mLightPotDistance;
extern bool CurrentLampType;

extern int shpID;
extern int bldNum;
extern std::vector<CultNode> totalShp;
extern std::vector<LocateBuilding> totalUserBld;
extern std::vector<HighCultureInfo> allHCI;
extern int RefreshTheConfigFile();
extern int SetGlobalVars();
extern void UpdateDialog01();
extern int MessageTo;
extern void MessageToDialog(char *msg);

tagLatLon LatLon;

DWORD WINAPI ThreadBatchCC (PVOID pParam);
HANDLE hLLThread;
int InfoLLFlag;

DWORD WINAPI ThreadBatchCC (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!BeginCreateCulture())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}


// CLatLon �Ի���

IMPLEMENT_DYNAMIC(CLatLon, CBkDialogST)

CLatLon::CLatLon(CWnd* pParent /*=NULL*/)
	: CBkDialogST(CLatLon::IDD, pParent)
	, m_listbox1(_T(""))
	, roadWide(20.0)
	, lightDistance(15.0)
	, woodCoverage(4000)
{

}

CLatLon::~CLatLon()
{
}

void CLatLon::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Control(pDX, 	IDC_TREE_CULTURE, m_tree);
	DDX_Control(pDX, 	IDC_LIST_CULTURE, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST_CULTURE, m_listbox1);
	DDX_Text(pDX, IDC_ROADWIDE, roadWide);
	DDV_MinMaxFloat(pDX, roadWide, 5, 1000);
	DDX_Text(pDX, IDC_LIGHTDISTANCE, lightDistance);
	DDV_MinMaxFloat(pDX, lightDistance, 0, 100);
	DDX_Text(pDX, IDC_WOODCOVERAGE, woodCoverage);
	DDV_MinMaxFloat(pDX, woodCoverage, 400, 2000000);
	DDX_Check(pDX, IDC_CHECK_LAMP, (int &)m_lampType);
}


BEGIN_MESSAGE_MAP(CLatLon, CBkDialogST)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE_CULTURE, &CLatLon::OnTvnSelchangedTree)
	ON_NOTIFY(TVN_ENDLABELEDIT, IDC_TREE_CULTURE, &CLatLon::OnTvnEndlabeleditTree)
	ON_NOTIFY(NM_DBLCLK, IDC_TREE_CULTURE, &CLatLon::OnDblClkTree)
	ON_BN_CLICKED(IDC_REFRESH_CULTURE, 	&CLatLon::OnBnClickedRefresh_Button)
	ON_BN_CLICKED(IDC_RESET_CULTURE, 		&CLatLon::OnBnClickedReset_Button)
	ON_BN_CLICKED(IDC_CREATE_BUTTON, &CLatLon::OnBnClickedCreateButton)
	ON_BN_CLICKED(IDC_CHECK_LAMP, &CLatLon::OnBnClickedLamp)

	ON_MESSAGE	 (WM_LATLON,   &CLatLon::OnListSet)
	ON_MESSAGE	 (WM_LATLON_UPDATE,  &CLatLon::OnUpdateEdit)

	ON_WM_CTLCOLOR()
	ON_WM_PAINT()

END_MESSAGE_MAP()


// CLatLon ��Ϣ�������
void CLatLon::OnBnClickedRefresh_Button()
{
	/* ��ʾ��Ϣ��*/
	MessageTo = 1;
	MessageToDialog("���� \"�Ļ���Ϣ\" ҳ��� \"ˢ��\" ��ť��");

	m_tree.DeleteAllItems();
	RefreshTheConfigFile();
	initTree();
}

void CLatLon::OnBnClickedReset_Button()
{
	/* ��ʾ��Ϣ��*/
	MessageTo = 1;
	MessageToDialog("���� \"�Ļ���Ϣ\" ҳ��� \"����\" ��ť��");

	shpID = 0;	totalShp.clear();
	bldNum = 0;	totalUserBld.clear();
	allHCI.clear();
	
	OnBnClickedRefresh_Button();
}




BOOL CLatLon::OnInitDialog()
{
	CDialog::OnInitDialog();

	typedef BOOL (WINAPI *lpfnSetLayeredWindowAttributes)(HWND hWnd, COLORREF crKey, BYTE bAlpha, DWORD dwFlags); 
	lpfnSetLayeredWindowAttributes SetLayeredWindowAttributes; 

	m_Refresh_Button.SubclassDlgItem(IDC_REFRESH_CULTURE, this);
	m_Refresh_Button.SetIcon(IDI_Go);
	m_Refresh_Button.DrawTransparent();

	m_Reset_Button.SubclassDlgItem(IDC_RESET_CULTURE, this);
	m_Reset_Button.SetIcon(IDI_Go);
	m_Reset_Button.DrawTransparent();
	
	m_Create_Button.SubclassDlgItem(IDC_CREATE_BUTTON, this);
	m_Create_Button.SetIcon(IDI_Go);
	m_Create_Button.DrawTransparent();

	initTree();

	m_hAccel=::LoadAccelerators(::AfxGetResourceHandle(),MAKEINTRESOURCE(IDR_ACCELERATOR1));

	/* ����һ�δ��벻֪����ʲô�õġ�*/
	COLORREF maskColor=RGB(0,0,0);

#if 1	//_DEBUG
	HMODULE hUser32 = GetModuleHandle("user32.dll"); //���ض�̬���ӿ�
#else
	HMODULE hUser32 = GetModuleHandle(L"user32.dll"); //���ض�̬���ӿ�
#endif
	SetLayeredWindowAttributes = (lpfnSetLayeredWindowAttributes)GetProcAddress(hUser32,"SetLayeredWindowAttributes"); 

	//ȡ��SetLayeredWindowAttributes����ָ�� 
	//Ϊ���ڼ���WS_EX_LAYERED��չ����
	SetWindowLong(this->GetSafeHwnd(), GWL_EXSTYLE, GetWindowLong(GetSafeHwnd(), GWL_EXSTYLE)^WS_EX_LAYERED); 

	//����SetLayeredWinowAttributes����
	SetLayeredWindowAttributes(this->GetSafeHwnd(), maskColor, 192, 2); 

	FreeLibrary(hUser32);   //�ͷŶ�̬���ӿ�

	// TODO:  �ڴ���Ӷ���ĳ�ʼ��
	UpdateData(false);
	// TODO:  �ڴ���Ӷ���ĳ�ʼ��

	return TRUE;  // return TRUE unless you set the focus to a control
	// �쳣: OCX ����ҳӦ���� FALSE
}

LRESULT CLatLon::OnListSet(WPARAM wParam, LPARAM lParam)
{
	//((CComboBox*)GetDlgItem(IDC_LIST1))->AddString(m_sCloudType[i]);
    char* temp=(char*)lParam;
	CString sTemp(temp);

	// Select the next item of the currently selected one.
	m_listbox1_ctrl.AddString(sTemp);
	int nCount = m_listbox1_ctrl.GetCount();
	int curSel = m_listbox1_ctrl.SetCurSel(nCount - 1);
	
	UpdateData(TRUE);
	return 0;
}

HBRUSH CLatLon::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  �ڴ˸��� DC ���κ�����
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
/*	case CTLCOLOR_EDIT:
		//	pDC ->SetBkColor(RGB(0,0,0));
		pDC->SetTextColor(RGB(0,0,0)); 
		hbr   =   ::CreateSolidBrush(RGB(188,197,230));   
		//pDC-> SetBkMode(TRANSPARENT); 
	return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); */
	}

	// TODO:  ���Ĭ�ϵĲ������軭�ʣ��򷵻���һ������
	return hbr;
}

void CLatLon::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: �ڴ˴������Ϣ����������
	// ��Ϊ��ͼ��Ϣ���� CDialog::OnPaint()
	CBitmap m_bmpBK;
	m_bmpBK.LoadBitmap(IDB_BITMAP1);

	CRect rect,rc;
	GetClientRect(&rect);//���Ŀ��ߴ磬�����ڿͻ���������
	GetWindowRect(&rc);

	//pParent-> ScreenToClient(  & rect   );

	BITMAP bitMap;//λͼ�ṹ��
	m_bmpBK.GetBitmap(&bitMap);//���ԭͼƬ�ߴ�

	CDC dcMem; //Ŀ��DC
	dcMem.CreateCompatibleDC(&dc); //������dc���ݵ��ڴ�DC
	dcMem.SelectObject(&m_bmpBK);//��λͼ����m_bmpBKѡ���ڴ�DC
	dc.StretchBlt(0,0,rect.Width(),rect.Height(),&dcMem,192,111,rect.Width(),rect.Height(),SRCCOPY);//bitMap.bmWidth-200,bitMap.bmHeight-250,SRCCOPY);

}

void CLatLon::OnBnClickedCreateButton()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData();

	/* ��ʾ��Ϣ��*/
	MessageTo = 1;
	MessageToDialog("���� \"�Ļ���Ϣ\" ҳ��� \"����\" ��ť��");

	CurrentLampType					 = m_lampType;		
	mRoadWidthInShp				   = roadWide;
	mLightPotDistance			   = lightDistance;
	mWoodCoverage	  			   = woodCoverage;

	GetDlgItem(IDC_CREATE_BUTTON)->EnableWindow(0);
	InfoLLFlag = 0;
	hLLThread = CreateThread(NULL, 0, ThreadBatchCC,(void*)&InfoLLFlag,0,NULL);
	GetDlgItem(IDC_CREATE_BUTTON)->EnableWindow(1);

}


void CLatLon::OnBnClickedLamp()
{
	if(((CButton*)GetDlgItem(IDC_CHECK_LAMP))-> GetCheck())
	{
		GetDlgItem(IDC_LIGHTDISTANCE)->EnableWindow(TRUE);
		GetDlgItem(IDC_STATIC13)->EnableWindow(TRUE);
	}
	else
	{
		GetDlgItem(IDC_LIGHTDISTANCE)->EnableWindow(FALSE);
		GetDlgItem(IDC_STATIC13)->EnableWindow(FALSE);
	}
}

LRESULT CLatLon::OnUpdateEdit(WPARAM wParam, LPARAM lParam)
{
	UpdateParameters();
	return 0;
}


void CLatLon::UpdateParameters()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(true);
}


void CLatLon::initTree()
{
	HTREEITEM hroot,hchild;
	char buf[256];
	int i,j;

	//���ڵ�����
	for(i=0;i<3;i++)
	{
		switch(i)
		{
		case 0:
			{
				sprintf(buf, "�û��Ļ���Ϣԭʼ�ļ�SHP (%d)", shpID);
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<shpID;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalShp[j].fileName, totalShp[j].pathName,totalShp[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//�ӽڵ�����
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 1:
			{
				sprintf(buf, "�û�������λ���FLT (%d)", bldNum);
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<bldNum;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalUserBld[j].fltFname, totalUserBld[j].fltPname,totalUserBld[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//�ӽڵ�����
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 2:
			{
				sprintf(buf, "�û��Ļ���Ϣ������ (%d)", allHCI.size());
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<allHCI.size();j++)
				{
					sprintf(buf, "%s    %9d   %2d", allHCI[j].texFname, allHCI[j].tileID, allHCI[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//�ӽڵ�����
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		default:
			break;
		}
	}

	m_tree.SetBkColor(RGB(247,247,255));
	m_tree.SetTextColor(RGB(0,0,255));
}


void CLatLon::OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	
	CString str=m_tree.GetItemText(pNMTreeView->itemNew.hItem);

	*pResult = 0;
}

void CLatLon::OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTVDISPINFO pTVDispInfo = reinterpret_cast<LPNMTVDISPINFO>(pNMHDR);
	
	HTREEITEM selItem,parentItem;
	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
		return;

	//��ȡ��������
	int indexTree,subIndex;
	subIndex = m_tree.GetItemData(selItem);
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	CString strText;
	m_tree.GetEditControl()->GetWindowText(strText.GetBuffer(200),200);
	
	char *strSel = (LPSTR)(LPCTSTR)strText;
	int strLen =  strlen(strSel);
	if(strLen > 0)
	{
		if((strSel[strLen-1] != '0') && (strSel[strLen-1] != '1'))
		{
			MessageBox("������ģ�������0����1");
			return;
		}
	}
	else
	{
		MessageBox("���ݴ���");
		return;
	}

	char *strItem = (LPSTR)(LPCTSTR)m_strItem;
	int lenItem = strlen(strItem);
	strItem[lenItem-1] = strSel[strLen-1];
	strText.Format("%s", strItem);
	m_tree.SetItemText(m_tree.GetSelectedItem(),strText);
	*pResult = 0;
}

void CLatLon::OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	HTREEITEM selItem,parentItem;

	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
	{
		return;
	}

	//��ȡ���ڵ�����
	int indexTree;
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	m_strItem = m_tree.GetItemText(selItem);
	m_tree.ModifyStyle(NULL,TVS_EDITLABELS);
	m_tree.EditLabel(m_tree.GetSelectedItem());
}

