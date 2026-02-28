/*  File: pepdisp.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * $Id: pepdisp.c,v 1.5 2017/06/05 20:59:30 mieg Exp $
 * Description: peptide display - using colcontrol.c
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 14 13:29 1998 (fw)
 * Created: Wed Jun  1 15:05:00 1994 (rd)
 *-------------------------------------------------------------------
 */

/*
#define PEP_DEBUG
*/

#include "acedb.h"

#include "freeout.h"
#include "display.h"
#include "lex.h"
#include "classes.h"
#include "sysclass.h"
#include "systags.h"
#include "tags.h"
#include "bs.h"
#include "peptide.h"
#include "a.h"
#include "pepdisp.h"
#include "dotter.h"
#include "graphAcedbInterface.h"

static MAPCONTROL pepMapDoMap(COLCONTROL control, KEY key, KEY from, int nmaps);

static Array allWindows = 0;

void recalcChecksums();

static void pepMapControlRemove(void *cp)
{ int i;
  COLCONTROL c = *((COLCONTROL *)cp);
  for ( i=0; i<arrayMax(allWindows); i++)
    if (arr(allWindows, i, COLCONTROL) == c)
      arr(allWindows, i, COLCONTROL) = 0;
}

/***************** configuration data ******************/
static void setActiveForKey(KEY key)
{
  MAPCONTROL map = currentMapControl();
  PEPLOOK *look = map->look;  

  look->block = FALSE;
  controlDraw();
}

static void toggleHeader()
{
  MAPCONTROL map = currentMapControl();
  COLCONTROL control = map->control;

  control->hideHeader = !control->hideHeader;
  controlDraw();
}

static void activeZoneGet()
{
  MAPCONTROL map = currentMapControl();
  PEPLOOK *look = map->look;  
  look->block = TRUE;
  displayBlock(setActiveForKey,"Double click key to set active zone");
}

static void none()
{
}


static COLOUROPT analysisButton[] =
{
  { (ColouredButtonFunc)none, 0, BLACK, WHITE, "Analysis...", 0 }
};

static struct ProtoStruct resSetProt ;


static struct ProtoStruct pepSequenceColumn = {
  0, pepSequenceCreate, 0, "pepSequence", 0, FALSE, 0, 0, 0
} ;

static struct ProtoStruct pepFeatureColumn = {
  0, pepFeatureCreate, 0, "pepFeature", 0, FALSE, FeatureConvert, 0, 0
} ;

static struct ProtoStruct hydrophobColumn = {
  0, hydrophobCreate, 0, "Hydrophobicity", 0, FALSE, 0, 0, 0
} ;

static struct ProtoStruct activeZoneColumn = {
  0, activeZoneCreate, 0, "pepActiveZone", 0, FALSE, 0, 0, 0
} ;

static struct ProtoStruct homolColumn = {
  0, homolCreate, 0, "Homol", 0, FALSE, homolConvert, 0, 0
} ;

static struct ProtoStruct homolNameColumn = {
  0, homolNameCreate, 0, "Homol_Name", 0, FALSE, homolConvert, 0, 0
} ;

static COLDEFAULT defaultColSet =
{
  { &mapLocatorColumn, TRUE,  "Locator" },
  { &mapScaleColumn, TRUE,    "Scale" },
  { &activeZoneColumn, TRUE,  "Active Zone" },
  { &resSetProt, FALSE,       "Stupid test" },
  { &pepSequenceColumn, TRUE, "Peptide sequence" },
  { &hydrophobColumn, TRUE,   "Hydrophobicity" },
  { &homolColumn, TRUE,       "Blastp" },
  { &homolColumn, TRUE,       "Blastx" },
  { &homolNameColumn, TRUE,   "Blastx names" },
  { &homolColumn, TRUE,       "Pfam-hmmls+Pfam-hmmfs" }, /* Don't use underscore since it's illegal to use in method names (esr) */
  { &homolColumn, FALSE,      "Pfam-hmmfs" },
  { 0, 0, 0 }
};

static COLPROTO pepProts[] = { 
  &mapLocatorColumn,
  &mapScaleColumn,
  &activeZoneColumn,
  &pepFeatureColumn,
  &pepSequenceColumn,
  &homolColumn,
  &homolNameColumn,
  &hydrophobColumn,
  0
} ;

static MENUOPT mainMenu[] =
{
  { graphDestroy,    "Quit" },
  { help,            "Help" },
  { graphPrint,      "Print Screen" },
  { displayPreserve, "Preserve" },
  { toggleHeader,    "Toggle Header" },
  { 0, 0 }
} ;

/************* little callback functions ***************/

static void pepFinalise (void *m)
{
  PEPLOOK *look = (PEPLOOK *)(m) ;
				/* what I do when this dies */
  if (arrayExists(look->pep))
    arrayDestroy (look->pep) ;
}
static void callDotter(void);

static void dumpActiveZone()
{
  MAPCONTROL map = currentMapControl();
  PEPLOOK *look = map->look;
  
  if(look->pep){
    /*BOOL result = */pepDumpFastA(look->pep,look->activeStart,look->activeEnd,name(map->key),0,0);
  }
}
static void setActiveZone(char *cp)
{
  MAPCONTROL map = currentMapControl();
  PEPLOOK *look = map->look;
  int x1,x2;
  BOOL changed = FALSE;

  if(sscanf(cp,"%d %d",&x1,&x2)==2){ /*x1-- ; x2-- ;  No zero mhmp 23.10.97*/
    if(x1 >= 0 && x1 <= arrayMax(look->pep) && x2 >= 0 && x2 <= arrayMax(look->pep)){
      if(x1 < x2){
	look->activeStart = x1;
	look->activeEnd = x2;
      }
      else{
	look->activeStart = x2;
	look->activeEnd = x1;
      }
      controlDraw();
      changed = TRUE;
    }
    else
      graphMessage(messprintf("integers must be between 0 and %d. You entered \"%s\"",arrayMax(look->pep),cp));
  }
  else{
    graphMessage(messprintf("format needs to be int int. Not \"%s\"",cp));
  }
  if(!changed){
    strncpy(look->activeText,messprintf("%d %d",look->activeStart,look->activeEnd),MAXACTIVETEXT);
    graphBoxDraw(look->activeregionBox,-1,-1);   
  }
}

static float header (MAPCONTROL map) 
{ 
  COLCONTROL control = map->control;
  PEPLOOK *look = map->look;
  static MENUOPT activeMenu[]=
    {
      { dumpActiveZone, "Fasta dump of Active zone" },
      { 0,0 }
    };
  float y;
  int x,box;

  look->aboutBox = graphBoxStart();
  graphTextFormat(BOLD);
  if (look->titleKey)
    { x = strlen (name(look->titleKey)) ;
      graphText (name(look->titleKey),1,0.5) ;
    }
  else
    { x = strlen (name(map->key)) ;
      graphText (name(map->key),1,0.5) ;
    }
  graphTextFormat (PLAIN_FORMAT) ;
  graphBoxEnd();
  box = graphButton("Active Zone:",activeZoneGet,x+2,0.5);
  graphBoxMenu(box, activeMenu);

  strncpy(look->activeText,messprintf("%d %d",look->activeStart,look->activeEnd),MAXACTIVETEXT);
  look->activeregionBox = graphBoxStart();
  graphTextEntry(look->activeText,MAXACTIVETEXT-1,x+15.5,0.5,setActiveZone);
  graphBoxEnd();
  graphBoxDraw(look->activeregionBox,BLACK,YELLOW);


  if (map->noButtons) /* may suppress header buttons for WWW */
    y = 1.4;
  else
    y = controlDrawButtons (map, x+28, 0.5, control->graphWidth) ; 
  
  look->messageBox = graphBoxStart();
  graphTextPtr(look->messageText,1,y + 0.25, MAXMESSAGETEXT-1);
  graphBoxEnd();
  graphBoxDraw (look->messageBox, BLACK, PALECYAN);

  return y + 2.0 ;			/* return number of lines used */
}

static void pepMapHeaderPick(MAPCONTROL map, int box,double x,double y)
{
  PEPLOOK *look = map->look;

  if(box == look->aboutBox)
      display(map->key, 0, TREE);
}

/*************** entry point ****************/

BOOL pepDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  COLCONTROL control ;
  Array pep ;
  KEYSET maps = keySetCreate();
  KEY _VManyPep, mm;
  BOOL haveControl = FALSE ;
  int i;


#if !defined(PEP_DEBUG)
  display (key, from, TREE) ;
  return FALSE ;
#endif

  lexaddkey("ManyPep", &mm,_VMainClasses);
  _VManyPep = KEYKEY(mm);

  /*
    lexaddkey("View", &mm, _VMainClasses);
    _VView = KEYKEY(mm);
  */
  if (!_VProtein)
    { KEY key ;
      lexaddkey ("Protein", &key, _VMainClasses) ; 
      _VProtein = KEYKEY(key) ;
    }

  if(class(key) != _VManyPep){
    if (class(key) == _VPeptide)
      lexword2key (name(key), &key, _VProtein) ;
    
    { OBJ obj ;
      if (class(key) != _VProtein &&
	  (obj = bsCreate (key)))
	{ bsGetKey (obj, str2tag ("Corresponding_protein"), &key) ;
	  bsDestroy (obj) ;
	}
    }
    
    if (!(pep = peptideGet (key)) || arrayMax(pep) < 2)
      { if (pep) arrayDestroy (pep) ;
	display (key, from, TREE) ;
	return FALSE ;
      }
    keySet(maps,keySetMax(maps))=key;
    from = key;
  }
  else{
    KEY k;
    OBJ obj = bsCreate(key);
    if(obj)
      { if(bsGetKey(obj,str2tag("Protein"),&k))
	  do {
	    OBJ obj2 ;
	    if (class(k) != _VProtein &&
		(obj2 = bsCreate (k)))
	      { bsGetKey (obj2, str2tag ("Corresponding_protein"), &k) ;
		bsDestroy (obj2) ;
	      }
    
	    if (!(pep = peptideGet (k)) || arrayMax(pep) < 2)
	      { if (pep) arrayDestroy (pep) ;
		display (k, from, TREE) ;
		return FALSE ;
	      }
	    else
	      arrayDestroy(pep);
	    keySet(maps,keySetMax(maps)) = k;
	  } while(bsGetKey(obj,_bsDown,&k));
	bsDestroy(obj);
      }
  }

  if (keySetMax(maps) == 0) /* failed to find any maps */
    { if (from) 
        display (from, 0, TREE) ;
      keySetDestroy(maps);

      return FALSE ;
    }

    /* Remove any non-graphic maps from the list and display them as trees */
  keySetSort(maps); /* for KeySetFind */
  for (i=0; i<arrayMax(maps); i++)
    { KEY map = arr(maps, i, KEY);
      OBJ obj = bsCreate(map);
      if (obj && bsFindTag(obj, _Non_graphic))
        { keySetRemove(maps, map);
          display(map, from, TREE);
        }
      bsDestroy(obj);
    }
  if (keySetMax(maps) == 0)
    { keySetDestroy(maps);
      return FALSE;
    }
  if (!allWindows)
    allWindows = arrayCreate(10, COLCONTROL);
  if (isOldGraph)
    { 
/* Stage two, see if we have any window's open with precisely the right
   maps already displayed, this is inspired by ace3. Note that if we have
   views tied down to any maps, we don't re-use, we re-draw */
  
      for (i=0 ; i<arrayMax(allWindows); i++)
        { int j;
          COLCONTROL c = arr(allWindows, i, COLCONTROL);
          if (!c)
            continue; /* deleted one */
          if (arrayMax(c->maps) != keySetMax(maps))
            continue; /* different number, failed */
          for (j=0; j<arrayMax(c->maps); j++)
            { MAPCONTROL m = arr(c->maps, j, MAPCONTROL);
              if (!keySetFind(maps, m->key, 0))
                break;

            }
          if (j == arrayMax(c->maps))
            break; /* found the match */
        }
  
      if (i != arrayMax(allWindows)) /* found match */
        { /*float centre;*/
          control = arr(allWindows, i, COLCONTROL);
          graphActivate(control->graph);
	  /*          for (i=0; i<arrayMax(control->maps); i++)
            { MAPCONTROL map = arr(control->maps, i, MAPCONTROL);
              if (getCentre(from, map->key, &centre))
                { map->centre = centre;
                  mapControlCursorSet(map, centre);
                }
            }*/
          haveControl = TRUE ;
        }
    }

  if (!haveControl)
/* Stage three make a new control, in new or existing window, then display
   the maps in it */
    { char buff[1000]; /* for the name */
      COLCONTROL *cp;
      buff[0] = 0;
      for (i=0; i<arrayMax(maps); i++)
        { KEY mk = arr(maps, i, KEY);
          strcat(buff, messprintf("%s ", name(mk)));
        }

      control = colControlCreate(isOldGraph, buff, PEPMAP);
      array(allWindows, arrayMax(allWindows), COLCONTROL) = control;
      /* The next bit ensures that the pointer to the control in allWindows
         is removed when it goes away */
      cp = (COLCONTROL *) handleAlloc(pepMapControlRemove, 
                                      control->handle,
                                      sizeof(COLCONTROL));
      *cp = control;
      
      for (i=0; i<arrayMax(maps); i++)
        { KEY mk = arr(maps, i, KEY);
            pepMapDoMap(control, mk, from, arrayMax(maps));
	}
    }
  
  control->from = from;
  controlDraw();
  keySetDestroy(maps);
  return TRUE;
}


static MAPCONTROL pepMapDoMap(COLCONTROL control, KEY key, KEY from, int nmaps )
{
  PEPLOOK *look ;
  MAPCONTROL map = mapControlCreate(control,0) ;
  static MENUOPT Analysis[]=
    {
      { callDotter, "Self-dotplot in Dotter " },
      { dumpActiveZone, "Fasta dump of Active zone" },
      { 0,0 }
    };

  look = (PEPLOOK*) halloc (sizeof(PEPLOOK), map->handle) ;
  blockSetFinalise ((char*)look, pepFinalise) ;
  /* must give map handle, since per map */
  map->look = (void*) look ;
  look->map = map ;
  

  if(from)
    look->titleKey = from;
  
  if(class(key) == str2tag("Peptide"))
    lexword2key(name(key), &key, _VProtein);

  if (!(look->pep = peptideGet (key)) || arrayMax(look->pep) < 2)
    { if (look->pep) arrayDestroy (look->pep) ;
      display (key, from, TREE) ;
      return map ;
    }
  
  
  
  *look->activeText = 0;

  map->displayName = PEPMAP ;		/* unused, but perhaps useful? */
  map->name = name(key) ;
  map->key = key ;
  map->menu = mainMenu ;
  map->headerPick = pepMapHeaderPick;
  map->drawHeader = header ;	/* to draw header */
  map->convertTrans = 0 ;	/* called each time you Draw */
  map->keyboard = 0 ;		/* called on keystrokes in your map */
  map->topMargin = 2 ;		/* size of header */
  map->prototypes = pepProts ;
  map->hasCursor = TRUE ;
  map->hasProjectionLines = TRUE ;
  map->viewTag = str2tag("Pepmap");
  map->min = 1;
  map->max = arrayMax (look->pep) ;
  map->centre = map->max/2 ;
  look->activeStart = 0;
  look->activeEnd = map->max;

  controlAddMap (control, map);
  controlAddButton(map,analysisButton ,Analysis);
  controlConvertPerm(map);
  /* mhmp 24.10.97 */
  map->mag = (control->graphHeight-4)/ (map->max - 1) ;	/* approx Whole */
  controlMakeColumns (map, defaultColSet, 0, FALSE) ;

  control->from = from ;
  /*  map->mag = (control->graphHeight-4)/ (map->max - 1) ;	approx Whole */
  controlDraw () ;
  return map ;
}

/*********************************************************************/
/************************ homol columns ******************************/

typedef struct { 
  KEY seq ;		/* subject, i.e. target sequence for match */
  KEY meth ;		/* type of match, e.g. BLASTP - controls display etc. */
  float score ;
  int qstart ;		/* query start and end, i.e. in self */
  int qend ;
  int sstart ;		/* subject start and end */
  int send ;
} HOMOL ;



static void callDotter(void)
{
  MAPCONTROL currentMap = currentMapControl();
  PEPLOOK *look = (PEPLOOK*)currentMap->look;
  char *seq, *cp;
  int i;
  
  if (!look->pep) return;
  
  seq = messalloc(currentMap->max+1);
  cp = seq;
  for (i = 0; i < currentMap->max; i++) *cp++ = pepDecodeChar[(int)(arr(look->pep, i, char))];
  dotter ('P',
          0,
          currentMap->name,
          seq,
          0,
          currentMap->name,
          seq,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
          0,
	  0);
}


/*****************************************************/
/**************   Gif Stuff **************************/

typedef struct pepGifStruct { int magic ; KEY key ; int x1, x2 ; Array pep ; } PEPGIF ;

void* pepGifGet (void *opaqueLook)
{
  KEY key ;
  Array pep = 0 ;
  int x1 = 1, x2 = 0 ;
  char *word ;
  PEPGIF* look = (PEPGIF*) opaqueLook ;

  if (look) pepGifDestroy (look) ;

  if (!(word = freeword()))
    goto usage ;
  if (!lexword2key (word, &key, _VProtein))
    { freeOutf ("// gif pepget error : Protein %s not known\n", word) ;
      goto usage ;
    }
  x1 = 1 ; x2 = 0 ;		/* default for whole sequence */
  if (freecheck ("w"))
    { word = freeword() ;
      if (!strcmp (word, "-coords"))
	{ if (!freeint (&x1) || !freeint (&x2) || (x2 <= x1) || x1 < 1)
	    goto usage ;
	}
      else
	goto usage ;
    }

  look = (PEPGIF *)messalloc (sizeof (struct pepGifStruct)) ;
  look->magic = 12345 ;
  look->key = key ;
  look->x1 = x1 ;
  look->x2 = x2 ;
  look->pep = pep ;
  return look ;

 usage:
  freeOutf ("// gif spepget error: usage: pepget protein [-coords x1 x2]  with 0 < x1 < x2\n") ;
  return 0 ;
}

/**************************************************************/
 
void pepGifSeq (void *opaqueLook)
{
  PEPGIF* look = (PEPGIF *) opaqueLook ;
  int i, x1, x2, nn = 0 ;
  char *cp, *cq, *buf = 0 ;

  if (!look ||  look->magic != 12345)
    { freeOutf ("// gif pepget error : no active protein\n") ;
      goto usage ;
    }

  if (!look->pep && !(look->pep = peptideGet (look->key)))
  { 
    freeOutf ("// gif pepget error : Protein sequence %s not known\n", name(look->key)) ;
    goto usage ;
  }
 x1 = look->x1 ; x2 = look->x2 ; nn = arrayMax(look->pep) ;
 if (nn && x2 == 0) look->x2 = x2 = nn ;
 if (x1 > nn) x1 = look->x1 = nn ;
 if (x2 > nn) x2 = look->x2 = nn ;
 if (x1 >= x2)
  { 
    freeOutf ("// gif pepget error : The protein has length=%d, smaller than your x1=%d\n", nn, look->x1) ;
    return ;
  }
 cp = arrp(look->pep, x1 - 1, char) ; 
 cq = buf = messalloc (x2 - x1 + 1) ;
 for (i = x1 ; i <= x2 ; i++, cp++, cq++)
   *cq = pepDecodeChar[(int)*cp] ;
 freeOutf (">%s\n%s\n",name(look->key), buf) ;
 messfree (buf) ;
    
 return ;
 usage:
  freeOutf ("// gif pepseq error: usage: pepget then pepseq\n") ;
}

/**************************************************************/
 
void pepGifAlign (void *opaqueLook)
{
  PEPGIF* look = (PEPGIF *) opaqueLook ;
  Array aa = 0, mypep = 0 ; 
  int i, j, ii ;
  char *cp, *cq, *buf = 0, *word ;
  int a1, a2, x1, x2, nn, nn2 ;
  KEY prot, type ;
  FILE *fil = 0 ; int level = 0 ;
  OBJ obj = 0 ;
  BSunit *uu ;

  if (!look ||  look->magic != 12345)
    { freeOutf ("// gif pepalign error : no active protein\n") ;
      goto usage ;
    }

  if (!look->pep && !(look->pep = peptideGet (look->key)))
    { freeOutf ("// gif pepalign error : Protein sequence %s not known\n", name(look->key)) ;
    goto usage ;
    }
  
  if (freecheck ("w") && (word = freeword()) &&  !strcmp (word, "-file"))
    { if (!freecheck ("w") || !(fil = filopen (freeword(), 0, "w")))
      goto usage ;
    level = freeOutSetFile (fil) ;
    }

  x1 = look->x1 ; x2 = look->x2 ; nn = arrayMax(look->pep) ;
  if (nn && x2 == 0) look->x2 = x2 = nn ;
  if (x1 > nn) x1 = look->x1 = nn ;
  if (x2 > nn) x2 = look->x2 = nn ;
  
  aa = arrayCreate (240, BSunit) ;
  if ((obj = bsCreate (look->key)))
   {
     bsGetArray (obj, str2tag("Homol"), aa, 8) ;
     bsDestroy (obj) ;
   }
  
  if (!arrayMax (aa))
    {
      freeOutf ("// gif pepalign error : No homology in %s\n", name(look->key)) ;
      arrayDestroy (aa) ;
      return ;
    }
 
 freeOutf ("%s/%d-%d ", name(look->key), look->x1, look->x2) ;
  cp = arrp(look->pep, x1 - 1, char) ; 
  cq = buf = messalloc (nn + 1) ;
  memset (buf, (int)'.', nn) ;
  buf[nn] = 0 ;
  for (i = x1 ; i <= x2 ; i++, cp++, cq++)
    *cq = pepDecodeChar[(int)*cp] ;
  freeOutf ("%s %s\n",buf, name(look->key)) ;
  messfree (buf) ;

  for (ii = 0 ; ii < arrayMax(aa) ; ii+= 8) 
    {
      uu = arrp(aa, ii, BSunit) ;
      type = uu[0].k ;
      prot = uu[1].k ;
      /*
	meth = uu[2].k ;
       score = uu[3].f ; 
      */
      a1 = uu[4].i ;
      a2 = uu[5].i ;
      x1 = uu[6].i ;
      x2 = uu[7].i ;
      mypep = 0 ;

      if (type == str2tag("Pep_homol") &&
	  (mypep = peptideGet(prot)))
	{
	  freeOutf ("%s/%d-%d ", name(prot), a1, a2, x1, x2) ;
	  cp = arrp(mypep, x1 - 1, char) ; 
	  cq = buf = messalloc (nn + 1) ;
	  memset (buf, (int)'.', nn) ;
	  buf[nn] = 0 ;
	  nn2 = arrayMax(mypep) ; cq += a1 - 1 ;
	  for (i = x1, j = a1 ; i <= x2 && i > 0 && i < nn2 ; i++, j++ , cp++, cq++)
	    if (j >= 0 && j < nn) 
	      *cq = pepDecodeChar[(int)*cp] ;
	  freeOutf ("%s %s\n",buf, name(prot)) ;
	  messfree (buf) ;
	  arrayDestroy (mypep) ;
	}
    }
  arrayDestroy (aa) ;
  if (level) { freeOutClose (level) ; filclose (fil) ; }
  return ;
usage:
  if (level) { freeOutClose (level) ; filclose (fil) ; }
  freeOutf ("// gif pepalign error: usage: pepget then pepalign\n") ;
}

/**************************************************************/

void pepGifDestroy (void *opaqueLook)
{
  PEPGIF* look = (PEPGIF *) opaqueLook ;

  if (look && look->magic == 12345)
    { 
      if (arrayExists (look->pep))
      arrayDestroy (look->pep) ;
    }
} 

/**************************************************************/


/*****************************************************/
/**************** end of file *****************/

 







 
 
 
 
 
 
 
