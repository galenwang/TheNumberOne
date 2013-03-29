// OperEnv.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "OperEnv.h"
#include "ProgressDialog.h"
#include "interface.h"
#include "travelmap.h"


// COperEnv �Ի���
extern TerrainInput currTi;
extern int needUpdate;
extern char baseTextDir[_MAX_PATH];
extern char baseHgtDir[_MAX_PATH];
extern char fsData[_MAX_PATH];
extern char workDirectory[_MAX_PATH];

extern int hTxtNum;
extern int hHgtNum;
extern unsigned int numHTI;
extern std::vector<HighNode> highTxtnode;
extern std::vector<HighNode> highHgtnode;
extern std::vector<HighTerrainInfo> allHTI;
extern int RefreshTheConfigFile();
extern int SetGlobalVars();
extern void UpdateDialog01();
extern void DoInitializationCheck();

DWORD WINAPI ThreadBatchCT (PVOID pParam);
HANDLE hThread;
int InfoFlag;
extern int MessageTo;
extern void MessageToDialog(char *msg);
extern string getCurrLODString();
extern int getCurrLodInt(const char *lodstr);

DWORD WINAPI ThreadBatchCT (PVOID pParam)
{
	int *flag=(int*)pParam;

	if(!BatchCreateTerrain())
	{
		*flag=2;
		return 0;
	}
	*flag=1;
	return 1;
}

CCharHandle::CCharHandle()
{
	//pResult=NULL;
}
CCharHandle::~CCharHandle()
{
	vector<char*>::iterator iter;
	for(iter=resV.begin();iter!=resV.end();iter++)
		delete *iter;
}

IMPLEMENT_DYNAMIC(COperEnv, CBkDialogST)

COperEnv::COperEnv(CWnd* pParent /*=NULL*/)
	: CBkDialogST(COperEnv::IDD, pParent)
	, m_Lod(_T(""))
	, m_listbox1(_T(""))
{


}

COperEnv::~COperEnv()
{
}



void COperEnv::DoDataExchange(CDataExchange* pDX)
{
	CBkDialogST::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TREE_TERRAIN, m_tree);
	DDX_Text(pDX, IDC_EDIT_LOD, m_Lod);
 	DDX_Control(pDX, IDC_LIST_TERRAIN, m_listbox1_ctrl);
	DDX_LBString(pDX, IDC_LIST_TERRAIN, m_listbox1);
}


BEGIN_MESSAGE_MAP(COperEnv, CBkDialogST)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE_TERRAIN, &COperEnv::OnTvnSelchangedTree)
	ON_NOTIFY(TVN_ENDLABELEDIT, IDC_TREE_TERRAIN, &COperEnv::OnTvnEndlabeleditTree)
	ON_NOTIFY(NM_DBLCLK, IDC_TREE_TERRAIN, &COperEnv::OnDblClkTree)
	ON_BN_CLICKED(IDC_REFRESH_TERRAIN, 	&COperEnv::OnBnClickedRefresh_Button)
	ON_BN_CLICKED(IDC_RESET_TERRAIN, 		&COperEnv::OnBnClickedReset_Button)
	ON_BN_CLICKED(IDC_PiLiang_Button, &COperEnv::OnBnClickedPiLiang_Button)
	ON_MESSAGE	 (WM_OPERENV,  &COperEnv::OnListSet)
	ON_MESSAGE	 (WM_OPERENV_UPDATE,  &COperEnv::OnUpdateEdit)
	ON_WM_CTLCOLOR()
	ON_WM_PAINT()
	ON_WM_TIMER()
END_MESSAGE_MAP()


// COperEnv ��Ϣ�������
void COperEnv::OnBnClickedRefresh_Button()
{
	/* ��ʾ��Ϣ��*/
	MessageTo = 0;
	MessageToDialog("���� \"�ؾ�����\" ҳ��� \"ˢ��\" ��ť��");

	m_tree.DeleteAllItems();
	RefreshTheConfigFile();
	initTree();
}

void COperEnv::OnBnClickedReset_Button()
{
	/* ��ʾ��Ϣ��*/
	MessageTo = 0;
	MessageToDialog("���� \"�ؾ�����\" ҳ��� \"����\" ��ť��");

	hTxtNum = 0;	highTxtnode.clear();
	hHgtNum = 0;	highHgtnode.clear();
	numHTI = 0;		allHTI.clear();
	
	OnBnClickedRefresh_Button();
}


void COperEnv::OnBnClickedPiLiang_Button()
{
	UpdateParameters();
	
	/* ��ʾ��Ϣ��*/
	MessageTo = 0;
	MessageToDialog("���� \"�ؾ�����\" ҳ��� \"����\" ��ť��");
	
	GetDlgItem(IDC_PiLiang_Button)->EnableWindow(0);

	InfoFlag = 0;
	hThread = CreateThread(NULL, 0, ThreadBatchCT,(void*)&InfoFlag,0,NULL);

	GetDlgItem(IDC_PiLiang_Button)->EnableWindow(1);
}

BOOL COperEnv::OnInitDialog()
{
	CDialog::OnInitDialog();

	typedef BOOL (WINAPI *lpfnSetLayeredWindowAttributes)(HWND hWnd, COLORREF crKey, BYTE bAlpha, DWORD dwFlags); 
	lpfnSetLayeredWindowAttributes SetLayeredWindowAttributes; 

	m_Refresh_Button.SubclassDlgItem(IDC_REFRESH_TERRAIN, this);
	m_Refresh_Button.SetIcon(IDI_Go);
	m_Refresh_Button.DrawTransparent();

	m_Reset_Button.SubclassDlgItem(IDC_RESET_TERRAIN, this);
	m_Reset_Button.SetIcon(IDI_Go);
	m_Reset_Button.DrawTransparent();
	
	m_PiLiang_Button.SubclassDlgItem(IDC_PiLiang_Button, this);
	m_PiLiang_Button.SetIcon(IDI_Go);
	m_PiLiang_Button.DrawTransparent();

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
	m_Lod.SetString(getCurrLODString().c_str());
	UpdateData(false);

	// ��һ����ʱ��
	this->SetTimer(1,500,NULL);//

	return TRUE;  // return TRUE unless you set the focus to a control
}


LRESULT COperEnv::OnListSet(WPARAM wParam, LPARAM lParam)
{
	char* temp=(char*)lParam;
	CString sTemp(temp);

	// Select the next item of the currently selected one.
	m_listbox1_ctrl.AddString(sTemp);
	int nCount = m_listbox1_ctrl.GetCount();
	int curSel = m_listbox1_ctrl.SetCurSel(nCount - 1);
	
	UpdateData(TRUE);
	return 0;
}


HBRUSH COperEnv::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  �ڴ˸��� DC ���κ�����
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
	}
	return hbr;
}

void COperEnv::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	BOOL first = true;

	CBitmap m_bmpBK;
	m_bmpBK.LoadBitmap(IDB_BITMAP1);

	CRect rect,rc;
	GetClientRect(&rect);//���Ŀ��ߴ磬�����ڿͻ���������
	GetWindowRect(&rc);

	BITMAP bitMap;//λͼ�ṹ��
	m_bmpBK.GetBitmap(&bitMap);//���ԭͼƬ�ߴ�

	CDC dcMem; //Ŀ��DC
	dcMem.CreateCompatibleDC(&dc); //������dc���ݵ��ڴ�DC
	dcMem.SelectObject(&m_bmpBK);//��λͼ����m_bmpBKѡ���ڴ�DC
	dc.StretchBlt(0,0,rect.Width(),rect.Height(),&dcMem,192,111,rect.Width(),rect.Height(),SRCCOPY);//bitMap.bmWidth-200,bitMap.bmHeight-250,SRCCOPY);
}

LRESULT COperEnv::OnUpdateEdit(WPARAM wParam, LPARAM lParam)
{
	UpdateParameters();
	return 0;
}


void COperEnv::UpdateParameters()
{
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	UpdateData(true);
	currTi.lod = getCurrLodInt(m_Lod.GetString());

}


void COperEnv::initTree()
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
				sprintf(buf, "�û��߾�������ͼƬTIF (%d)", hTxtNum);
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hTxtNum;j++)
				{
					sprintf(buf, "%s    %s", highTxtnode[j].fileName,highTxtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//�ӽڵ�����
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 1:
			{
				sprintf(buf, "�û��߾��ȸ߳�ͼƬTIF (%d)", hHgtNum);
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hHgtNum;j++)
				{
					sprintf(buf, "%s    %s", highHgtnode[j].fileName, highHgtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//�ӽڵ�����
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 2:
			{
				sprintf(buf, "�û��߾��ȵؾ�������� (%d)", numHTI);
				hroot=m_tree.InsertItem(buf);
				//���ڵ�����
				m_tree.SetItemData(hroot,i);
				for(j=0;j<numHTI;j++)
				{
					sprintf(buf, "%s    %s    %d  %2d", allHTI[j].texFname, allHTI[j].hgtFname, allHTI[j].tileID, allHTI[j].isOK);
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


void COperEnv::OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	
	CString str=m_tree.GetItemText(pNMTreeView->itemNew.hItem);

	*pResult = 0;
}

void COperEnv::OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult)
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

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("�����ݲ�������ģ�");
		return;
	}

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

	if(indexTree == 2)
	{
		allHTI[subIndex].isOK = atoi(&strSel[strLen-1]);
	}

	char *strItem = (LPSTR)(LPCTSTR)m_strItem;
	int lenItem = strlen(strItem);
	strItem[lenItem-1] = strSel[strLen-1];
	strText.Format("%s", strItem);
	m_tree.SetItemText(m_tree.GetSelectedItem(),strText);
	*pResult = 0;
}

void COperEnv::OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult)
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

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("�����ݲ�������ģ�");
		return;
	}

	m_strItem = m_tree.GetItemText(selItem);
	m_tree.ModifyStyle(NULL,TVS_EDITLABELS);
	m_tree.EditLabel(m_tree.GetSelectedItem());
}

void COperEnv::OnTimer(UINT_PTR nIDEvent)
{
	// ִ�г�ʼ�����
	DoInitializationCheck();

	CDialog::OnTimer(nIDEvent);
	this->KillTimer(1); 
}

