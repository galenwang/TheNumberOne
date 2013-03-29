// ManageProfile.cpp : 实现文件
//

#include "stdafx.h"
#include "dijingzhizuo.h"
#include "ManageProfile.h"
#include "interface.h"
#include "travelmap.h"

extern int hTxtNum;
extern int hHgtNum;
extern unsigned int numHTI;
extern int cultID;
extern int shpID;
extern int bldNum;

extern std::vector<HighNode> highTxtnode;
extern std::vector<HighNode> highHgtnode;
extern std::vector<HighTerrainInfo> allHTI;
extern std::vector<CultNode> totalcult;
extern std::vector<CultNode> totalShp;
extern std::vector<LocateBuilding> totalUserBld;
extern int RefreshTheConfigFile();
extern int SetGlobalVars();
extern void UpdateDialog01();
// CManageProfile 对话框

IMPLEMENT_DYNAMIC(CManageProfile, CDialog)

CManageProfile::CManageProfile(CWnd* pParent /*=NULL*/)
	: CDialog(CManageProfile::IDD, pParent)
{
}

CManageProfile::~CManageProfile()
{
}

void CManageProfile::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TREE1, m_tree);
}


BEGIN_MESSAGE_MAP(CManageProfile, CDialog)
	ON_WM_CTLCOLOR()
	ON_WM_PAINT()
	ON_WM_ERASEBKGND()
	ON_BN_CLICKED(IDC_EDIT, &CManageProfile::OnBnClickedEdit)
	ON_NOTIFY(TVN_SELCHANGED, IDC_TREE1, &CManageProfile::OnTvnSelchangedTree)
	ON_NOTIFY(TVN_ENDLABELEDIT, IDC_TREE1, &CManageProfile::OnTvnEndlabeleditTree)
	ON_NOTIFY(NM_DBLCLK, IDC_TREE1, &CManageProfile::OnDblClkTree)
	ON_BN_CLICKED(IDC_SAVE, &CManageProfile::OnBnClickedSave)
	ON_BN_CLICKED(IDC_DELETE, &CManageProfile::OnBnClickedDelete)
	ON_BN_CLICKED(IDC_REFRESH, &CManageProfile::OnBnClickedRefresh)
END_MESSAGE_MAP()


// CManageProfile 消息处理程序
HBRUSH CManageProfile::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
	HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

	// TODO:  如果默认的不是所需画笔，则返回另一个画笔
	switch(nCtlColor)
	{
	case CTLCOLOR_STATIC:
		pDC ->SetBkMode(TRANSPARENT);
		pDC->SetTextColor(RGB(0,0,0)); 
		return (HBRUSH)::GetStockObject(HOLLOW_BRUSH); 
	}
	return hbr;
}

void CManageProfile::OnPaint()
{
	CPaintDC dc(this); // device context for painting
	// TODO: 在此处添加消息处理程序代码
	// 不为绘图消息调用 CDialog::OnPaint()
	CBitmap m_bmpBK;
	m_bmpBK.LoadBitmap(IDB_BITMAP1);

	CRect rect,rc;
	GetClientRect(&rect);//获得目标尺寸，即窗口客户区的坐标
	GetWindowRect(&rc);

	//pParent-> ScreenToClient(  & rect   );

	BITMAP bitMap;//位图结构体
	m_bmpBK.GetBitmap(&bitMap);//获得原图片尺寸

	CDC dcMem; //目标DC
	dcMem.CreateCompatibleDC(&dc); //创建与dc兼容的内存DC
	dcMem.SelectObject(&m_bmpBK);//将位图对象m_bmpBK选入内存DC
	dc.StretchBlt(0,0,rect.Width(),rect.Height(),&dcMem,192,111,rect.Width(),rect.Height(),SRCCOPY);//bitMap.bmWidth-200,bitMap.bmHeight-250,SRCCOPY);

}

BOOL CManageProfile::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	return true;//CDialog::OnEraseBkgnd(pDC);
}

BOOL CManageProfile::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  在此添加额外的初始化
	m_Edit_Button.SubclassDlgItem(IDC_PROFILE, this);
	m_Edit_Button.SetIcon(IDI_Go);
	m_Edit_Button.DrawTransparent();

	m_Save_Button.SubclassDlgItem(IDC_SAVE, this);
	m_Save_Button.SetIcon(IDI_Go);
	m_Save_Button.DrawTransparent();

	m_Delete_Button.SubclassDlgItem(IDC_DELETE, this);
	m_Delete_Button.SetIcon(IDI_Go);
	m_Delete_Button.DrawTransparent();

	m_Refresh_Button.SubclassDlgItem(IDC_REFRESH, this);
	m_Refresh_Button.SetIcon(IDI_Go);
	m_Refresh_Button.DrawTransparent();
	
	initTree();

	m_hAccel=::LoadAccelerators(::AfxGetResourceHandle(),MAKEINTRESOURCE(IDR_ACCELERATOR1));

	return TRUE;  // return TRUE unless you set the focus to a control
	// 异常: OCX 属性页应返回 FALSE
}

//树内容初始化
void CManageProfile::initTree()
{
	HTREEITEM hroot,hchild;
	char buf[256];
	int i,j;

	//父节点内容
	for(i=0;i<6;i++)
	{
		switch(i)
		{
		case 0:
			{
				sprintf(buf, "用户高精度纹理图片TIF (%d)", hTxtNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hTxtNum;j++)
				{
					sprintf(buf, "%s    %s", highTxtnode[j].fileName,highTxtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 1:
			{
				sprintf(buf, "用户高精度高程图片TIF (%d)", hHgtNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<hHgtNum;j++)
				{
					sprintf(buf, "%s    %s", highHgtnode[j].fileName, highHgtnode[j].pathName);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 2:
			{
				sprintf(buf, "用户高精度地景制作情况 (%d)", numHTI);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<numHTI;j++)
				{
					sprintf(buf, "%s    %s  %2d", allHTI[j].texFname, allHTI[j].hgtFname,allHTI[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 3:
			{
				sprintf(buf, "用户文化信息制作情况AC (%d)", cultID);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<cultID;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalcult[j].fileName, totalcult[j].pathName,totalcult[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 4:
			{
				sprintf(buf, "用户文化信息制作情况SHP (%d)", shpID);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<shpID;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalShp[j].fileName, totalShp[j].pathName,totalShp[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
					m_tree.SetItemData(hchild,j);
				}
			}
			break;
		case 5:
			{
				sprintf(buf, "用户机场定位情况FLT (%d)", bldNum);
				hroot=m_tree.InsertItem(buf);
				//父节点索引
				m_tree.SetItemData(hroot,i);
				for(j=0;j<bldNum;j++)
				{
					sprintf(buf, "%s    %s   %2d", totalUserBld[j].fltFname, totalUserBld[j].fltPname,totalUserBld[j].isOK);
					hchild=m_tree.InsertItem(buf,hroot);
					//子节点索引
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


void CManageProfile::OnBnClickedEdit()
{
	// TODO: 在此添加控件通知处理程序代码
	HTREEITEM selItem,parentItem;

	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
	{
		return;
	}

	//获取父节点索引
	int indexTree;
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("该内容不允许更改！");
		return;
	}

	m_strItem = m_tree.GetItemText(selItem);

	m_tree.ModifyStyle(NULL,TVS_EDITLABELS);
	m_tree.EditLabel(m_tree.GetSelectedItem());
}

void CManageProfile::OnTvnSelchangedTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTREEVIEW pNMTreeView = reinterpret_cast<LPNMTREEVIEW>(pNMHDR);
	// TODO: 在此添加控件通知处理程序代码
	CString str=m_tree.GetItemText(pNMTreeView->itemNew.hItem);

	*pResult = 0;
}

void CManageProfile::OnTvnEndlabeleditTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMTVDISPINFO pTVDispInfo = reinterpret_cast<LPNMTVDISPINFO>(pNMHDR);
	// TODO: 在此添加控件通知处理程序代码
	HTREEITEM selItem,parentItem;
	selItem = m_tree.GetSelectedItem();
	if(m_tree.ItemHasChildren(selItem))
		return;

	//获取父子索引
	int indexTree,subIndex;
	subIndex = m_tree.GetItemData(selItem);
	parentItem = m_tree.GetParentItem(selItem);
	if(parentItem == NULL) return;
	indexTree = m_tree.GetItemData(parentItem);

	if((indexTree == 0) || (indexTree == 1))
	{
		MessageBox("该内容不允许更改！");
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
			MessageBox("错误更改！请输入0或者1");
			return;
		}
	}
	else
	{
		MessageBox("内容错误！");
		return;
	}

	if(indexTree == 2)
	{
		allHTI[subIndex].isOK = atoi(&strSel[strLen-1]);
	}
	else if(indexTree == 3)
	{
		totalcult[subIndex].isOK = atoi(&strSel[strLen-1]);
	}
	else if(indexTree == 4)
	{
		totalShp[subIndex].isOK = atoi(&strSel[strLen-1]);
	}
	else if(indexTree == 5)
	{
		totalUserBld[subIndex].isOK = atoi(&strSel[strLen-1]);
	}

	char *strItem = (LPSTR)(LPCTSTR)m_strItem;
	int lenItem = strlen(strItem);

	strItem[lenItem-1] = strSel[strLen-1];

	strText.Format("%s", strItem);

	m_tree.SetItemText(m_tree.GetSelectedItem(),strText);

	*pResult = 0;
}

void CManageProfile::OnDblClkTree(NMHDR *pNMHDR, LRESULT *pResult)
{
	OnBnClickedEdit();
}

void CManageProfile::OnBnClickedSave()
{
	// TODO: 在此添加控件通知处理程序代码
	CheckConfig(1);
	m_tree.DeleteAllItems();
	initTree();
	MessageBox("所有更改已保存！");
}

void CManageProfile::OnBnClickedDelete()
{
	// TODO: 在此添加控件通知处理程序代码
	HTREEITEM selItem,parentItem;
	selItem = m_tree.GetSelectedItem();
	if(selItem)
	{
		//树节点索引
		int indexParent,indexChild,i;

		indexParent = -1;
		indexChild = -1;
		//子节点
		if(m_tree.ItemHasChildren(selItem))
		{
			indexParent = m_tree.GetItemData(selItem);
		}
		else
		{
			parentItem = m_tree.GetParentItem(selItem);
			if(parentItem == NULL) return;
			indexParent = m_tree.GetItemData(parentItem);
			//if(indexParent == NULL) return;
			indexChild = m_tree.GetItemData(selItem);
		}
	
		switch(indexParent)
		{
			case 0:
				{
					if(indexChild == -1)
					{
						hTxtNum = 0;	highTxtnode.clear();
					}
					else
					{						
						for(i=indexChild;i<(hTxtNum-1);i++)
						{
							highTxtnode[i] = highTxtnode[i+1];
						}
						hTxtNum--;	highTxtnode.pop_back();
					}				
				}
				break;
			case 1:
				{
					if(indexChild == -1)
					{
						hHgtNum = 0;	highHgtnode.clear();
					}
					else
					{						
						for(i=indexChild;i<(hHgtNum-1);i++)
						{
							highHgtnode[i] = highHgtnode[i+1];
						}
						hHgtNum--;	highHgtnode.pop_back();
					}					
				}
				break;
			case 2:
				{
					if(indexChild == -1)
					{
						numHTI = 0;		allHTI.clear();
					}
					else
					{						
						for(i=indexChild;i<(numHTI-1);i++)
						{
							allHTI[i] = allHTI[i+1];
						}
						numHTI--;	allHTI.pop_back();
					}	
				}
				break;
			case 3:
				{
					if(indexChild == -1)
					{
						cultID = 0;		totalcult.clear();
					}
					else
					{						
						for(i=indexChild;i<(cultID-1);i++)
						{
							totalcult[i] = totalcult[i+1];
						}
						cultID--;	totalcult.pop_back();
					}	
				}
				break;
			case 4:
				{
					if(indexChild == -1)
					{
						shpID = 0;		totalShp.clear();
					}
					else
					{						
						for(i=indexChild;i<(shpID-1);i++)
						{
							totalShp[i] = totalShp[i+1];
						}
						shpID--;	totalShp.pop_back();
					}	
				}
				break;
			case 5:
				{
					if(indexChild == -1)
					{
						bldNum = 0;		totalUserBld.clear();
					}
					else
					{						
						for(i=indexChild;i<(bldNum-1);i++)
						{
							totalUserBld[i] = totalUserBld[i+1];
						}
						bldNum--;	totalUserBld.pop_back();
					}	
				}
				break;
			default:
				break;
		}
		m_tree.DeleteItem(selItem);
	}
}

void CManageProfile::OnBnClickedRefresh()
{
	// TODO: 在此添加控件通知处理程序代码
	/* 将路径更新一下 */
	UpdateDialog01();
	SetGlobalVars();

	m_tree.DeleteAllItems();
	RefreshTheConfigFile();
	initTree();
}
