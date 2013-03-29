#ifndef _MAP_GLOBAL_H_
#define _MAP_GLOBAL_H_

#define infinity 1e20           //定义为无穷大
#define EP 1e-10                //定义无穷小   

typedef unsigned char uchar;
typedef unsigned int  uint;

typedef struct DrawParam {
	
	double m_StartX;       //当前坐标原点横坐标
	double m_StartY;       //当前坐标原点纵坐标 
	int m_ScreenWidth;  //当前屏幕宽度
	int m_ScreenHeigh;  //当前屏幕高度
	double m_Top;
	double m_left;
	float m_Scale;      //当前绘图比例尺     


} DRAWPARAM;

/************************************************************************/
/*                              SwapWord()                              */
/*                                                                      */
/*      Swap a 2, 4 or 8 byte word.                                     */
/************************************************************************/

 static void	SwapWord( int length, void * wordP )

{
    int		i;
    uchar	temp;

    for( i=0; i < length/2; i++ )
    {
	temp = ((uchar *) wordP)[i];
	((uchar *)wordP)[i] = ((uchar *) wordP)[length-i-1];
	((uchar *) wordP)[length-i-1] = temp;
    }
}

#endif //_MAP_GLOBAL_H_