
#ifndef __PL_DLL_H
#define __PL_DLL_H
#ifdef USETHISz
// For an unknown compiler or static built we clear the macros
#ifndef PLDLLEXPORT
  #define PLDLLEXPORT
  #define PLDLLIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case), if
// cmake didn't set the corresponding macro xxxx_EXPORTS when the
// corresponding library is built (DLLIMPEXP is set to DLLEXPORT
// then)
#if defined ( plplotd_EXPORTS )
  #define PLDLLIMPEXP    PLDLLEXPORT
  #define PLDLLIMPEXP_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP    PLDLLIMPORT
  #define PLDLLIMPEXP_DATA( type )    PLDLLIMPORT type
#endif

// for dynamic drivers set the macros correctly. If a shared library is built,
// but dyanmic drivers disabled, the driver dll macros are the same as the#endif
// plplot dll macros

  #define PLDLLIMPEXP_DRIVER    PLDLLIMPEXP
  #define PLDLLIMPEXP_DRIVER_DATA( type )    PLDLLIMPEXP_DATA( type )


#if defined ( plplotcxxd_EXPORTS )
  #define PLDLLIMPEXP_CXX    PLDLLEXPORT
  #define PLDLLIMPEXP_CXX_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_CXX    PLDLLIMPORT
  #define PLDLLIMPEXP_CXX_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplotf95cd_EXPORTS )
  #define PLDLLIMPEXP_F95C    PLDLLEXPORT
  #define PLDLLIMPEXP_F95C_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_F95C    PLDLLIMPORT
  #define PLDLLIMPEXP_F95C_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplotwxwidgetsd_EXPORTS )
  #define PLDLLIMPEXP_WX    PLDLLEXPORT
  #define PLDLLIMPEXP_WX_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_WX    PLDLLIMPORT
  #define PLDLLIMPEXP_WX_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( tclmatrixd_EXPORTS )
  #define PLDLLIMPEXP_TCLMAT    PLDLLEXPORT
  #define PLDLLIMPEXP_TCLMAT_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_TCLMAT    PLDLLIMPORT
  #define PLDLLIMPEXP_TCLMAT_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplottcltk_Maind_EXPORTS ) | defined ( plplottcltkd_EXPORTS )
  #define PLDLLIMPEXP_TCLTK    PLDLLEXPORT
  #define PLDLLIMPEXP_TCLTK_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_TCLTK    PLDLLIMPORT
  #define PLDLLIMPEXP_TCLTK_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplotgnome2d_EXPORTS )
  #define    PLDLLEXPORT
  #define PLDLLIMPEXP_GNOME2_DATA( type )    PLDLLEXPORT type
#else
  #define    PLDLLIMPORT
  #define PLDLLIMPEXP_GNOME2_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( cplplotcanvasmodule_EXPORTS )
  #define PLDLLIMPEXP_CPLPLOTCANVASMODULE    PLDLLEXPORT
  #define PLDLLIMPEXP_CPLPLOTCANVASMODULE_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_CPLPLOTCANVASMODULE    PLDLLIMPORT
  #define PLDLLIMPEXP_CPLPLOTCANVASMODULE_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplot_widgetmodule_EXPORTS )
  #define PLDLLIMPEXP_PLPLOT_WIDGETMODULE    PLDLLEXPORT
  #define PLDLLIMPEXP_PLPLOT_MODULE_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_PLPLOT_MODULE          PLDLLIMPORT
  #define PLDLLIMPEXP_PLPLOT_MODULE_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplotqtd_EXPORTS )
  #define PLDLLIMPEXP_QT    PLDLLEXPORT
  #define PLDLLIMPEXP_QT_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_QT    PLDLLIMPORT
  #define PLDLLIMPEXP_QT_DATA( type )    PLDLLIMPORT type
#endif

#if defined ( plplot_pyqt4_EXPORTS )
  #define PLDLLIMPEXP_PYQT4    PLDLLEXPORT
  #define PLDLLIMPEXP_PYQT4_DATA( type )    PLDLLEXPORT type
#else
  #define PLDLLIMPEXP_PYQT4    PLDLLIMPORT
  #define PLDLLIMPEXP_PYQT4_DATA( type )    PLDLLIMPORT type
#endif
#endif
#endif // __PL_DLL_H
