#ifndef BigModelDebug_h___
#define BigModelDebug_h___


#ifdef DEBUG_ON
	static unsigned __active_debug_flags__ = 0;
#endif


#define DEBUG_INPUT_DATA		0x00000001
#define DEBUG_YEARLY_DATA		0x00000002
#define DEBUG_FIRE_DATA			0x00000004


#ifdef DEBUG_ON
#define DEBUG(flag,output) \
{ if ( flag & __active_debug_flags__ ) \
	cout << output << endl; }
#else
#define DEBUG(flag,output) ;
#endif

#endif




























