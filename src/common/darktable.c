/*
    This file is part of darktable,
    copyright (c) 2009--2011 johannes hanika.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "common/darktable.h"
#include "common/collection.h"
#include "common/exif.h"
#include "common/fswatch.h"
#include "common/grealpath.h"
#include "common/pwstorage/pwstorage.h"
#ifdef HAVE_GPHOTO2
#include "common/camera_control.h"
#endif
#include "common/film.h"
#include "common/image.h"
#include "common/image_cache.h"
#include "common/imageio_module.h"
#include "common/points.h"
#include "common/opencl.h"
#include "common/utility.h"
#include "develop/imageop.h"
#include "develop/blend.h"
#include "libs/lib.h"
#include "views/view.h"
#include "control/control.h"
#include "control/conf.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

#if !defined(__APPLE__) && !defined(__FreeBSD__) || defined(__WIN32__)
#include <malloc.h>
#endif
#ifdef __APPLE__
#include <sys/malloc.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

darktable_t darktable;
const char dt_supported_extensions[] = "3fr,arw,bay,bmq,cap,cine,cr2,crw,cs1,dc2,dcr,dng,erf,fff,exr,ia,iiq,jpg,jpeg,k25,kc2,kdc,mdc,mef,mos,mrw,nef,nrw,orf,pef,pfm,pxn,qtk,raf,raw,rdc,rw2,rwl,sr2,srf,sti,tif,tiff,x3f";

static int usage(const char *argv0)
{
  printf("usage: %s [-d {all,cache,control,dev,fswatch,camctl,perf,pwstorage,opencl,sql}] [IMG_1234.{RAW,..}|image_folder/]", argv0);
#ifdef HAVE_OPENCL
  printf(" [--disable-opencl]");
#endif
  printf("\n");
  return 1;
}

typedef void (dt_signal_handler_t)(int) ;
static dt_signal_handler_t *_dt_sigill_old_handler = NULL;

static
void _dt_sigill_handler(int param)
{
  GtkWidget *dlg = gtk_message_dialog_new(NULL,GTK_DIALOG_MODAL,GTK_MESSAGE_ERROR,GTK_BUTTONS_OK,_(
"darktable has trapped an illegal instruction which probably means that \
an invalid processor optimized codepath is used for your cpu, please try reproduce the crash running 'gdb darktable' from \
the console and post the backtrace log to mailing list with information about your CPU and where you got the package from."));
  gtk_dialog_run(GTK_DIALOG(dlg));
  
  /* pass it further the old handler*/
  _dt_sigill_old_handler(param);
}

#if defined(__i386__) && defined(__PIC__)
#define cpuid(level, a, b, c, d) \
  __asm__ ("xchgl %%ebx, %1\n" \
        "cpuid\n" \
        "xchgl  %%ebx, %1\n" \
        : "=a" (a), "=r" (b), "=c" (c), "=d" (d)	\
        : "0" (level) \
      )
#else
#define cpuid(level, a, b, c, d) \
  __asm__ ("cpuid"	\
    : "=a" (a), "=b" (b), "=c" (c), "=d" (d) \
    : "0" (level) \
    )
#endif

static 
void dt_check_cpu(int argc,char **argv) 
{
  /* hook up SIGILL handler */
  _dt_sigill_old_handler = signal(SIGILL,&_dt_sigill_handler);

  /* call cpuid for  SSE level */
  int ax,bx,cx,dx;
  cpuid(0x1,ax,bx,cx,dx);
  
  ax = bx = 0;
  char message[512]={0};
  strcat(message,_("SIMD extensions found: "));
  if((cx & 1) && (darktable.cpu_flags |= DT_CPU_FLAG_SSE3))
    strcat(message,"SSE3 ");
  if( ((dx >> 26) & 1) && (darktable.cpu_flags |= DT_CPU_FLAG_SSE2))
    strcat(message,"SSE2 ");
  if (((dx >> 25) & 1) && (darktable.cpu_flags |= DT_CPU_FLAG_SSE))
   strcat(message,"SSE ");
  if (!darktable.cpu_flags)
    strcat(message,"none");
 
  /* for now, bail out if SSE2 is not availble */
  if(!(darktable.cpu_flags & DT_CPU_FLAG_SSE2))
  {
    gtk_init (&argc, &argv);
    
    GtkWidget *dlg = gtk_message_dialog_new(NULL,GTK_DIALOG_MODAL,GTK_MESSAGE_ERROR,GTK_BUTTONS_OK,_(
"darktable is very cpu intensive and uses SSE2 SIMD instructions \
for heavy calculations. This gives a better user experience but also defines a minimum \
processor requirement.\n\nThe processor in YOUR system does NOT support SSE2. \
darktable will now close down.\n\n%s"),message);
    
    gtk_dialog_run(GTK_DIALOG(dlg));
    
    exit(11);
  }
}

static void strip_semicolons_from_keymap(const char* path)
{
  char pathtmp[1024];
  FILE *fin = fopen(path, "r");
  FILE *fout;
  int i;
  int c = '\0';

  snprintf(pathtmp, 1024, "%s_tmp", path);
  fout = fopen(pathtmp, "w");

  // First ignoring the first three lines
  for(i = 0; i < 3; i++)
  {
    c = fgetc(fin);
    while(c != '\n')
      c = fgetc(fin);
  }

  // Then ignore the first two characters of each line, copying the rest out
  while(c != EOF)
  {
    fseek(fin, 2, SEEK_CUR);
    do
    {
      c = fgetc(fin);
      if(c != EOF)
        fputc(c, fout);
    }while(c != '\n' && c != EOF);
  }

  fclose(fin);
  fclose(fout);
  g_file_delete(g_file_new_for_path(path), NULL, NULL);
  g_file_move(g_file_new_for_path(pathtmp), g_file_new_for_path(path), 0,
              NULL, NULL, NULL, NULL);
}

int dt_init(int argc, char *argv[], const int init_gui)
{
#ifndef __SSE2__
  fprintf(stderr, "[dt_init] unfortunately we depend on SSE2 instructions at this time.\n");
  fprintf(stderr, "[dt_init] please contribute a backport patch (or buy a newer processor).\n");
  return 1;
#endif
#ifdef M_MMAP_THRESHOLD
  mallopt(M_MMAP_THRESHOLD,128*1024) ; /* use mmap() for large allocations */   
#endif
  bindtextdomain (GETTEXT_PACKAGE, DARKTABLE_LOCALEDIR);
  bind_textdomain_codeset (GETTEXT_PACKAGE, "UTF-8");
  textdomain (GETTEXT_PACKAGE);

  
  // init all pointers to 0:
  memset(&darktable, 0, sizeof(darktable_t));

  darktable.progname = argv[0];
   
  // database
  gchar *dbfilenameFromCommand = NULL;

  darktable.num_openmp_threads = 1;
#ifdef _OPENMP
  darktable.num_openmp_threads = omp_get_num_procs();
#endif
  darktable.unmuted = 0;
  char *image_to_load = NULL;
  for(int k=1; k<argc; k++)
  {
    if(argv[k][0] == '-')
    {
      if(!strcmp(argv[k], "--help"))
      {
        return usage(argv[0]);
      }
      else if(!strcmp(argv[k], "--version"))
      {
        printf("this is "PACKAGE_STRING"\ncopyright (c) 2009-2011 johannes hanika\n"PACKAGE_BUGREPORT"\n");
        return 1;
      }
      else if(!strcmp(argv[k], "--library"))
      {
        dbfilenameFromCommand = argv[++k];
      }
      else if(argv[k][1] == 'd' && argc > k+1)
      {
        if(!strcmp(argv[k+1], "all"))            darktable.unmuted = 0xffffffff;   // enable all debug information
        else if(!strcmp(argv[k+1], "cache"))     darktable.unmuted |= DT_DEBUG_CACHE;   // enable debugging for lib/film/cache module
        else if(!strcmp(argv[k+1], "control"))   darktable.unmuted |= DT_DEBUG_CONTROL; // enable debugging for scheduler module
        else if(!strcmp(argv[k+1], "dev"))       darktable.unmuted |= DT_DEBUG_DEV; // develop module
        else if(!strcmp(argv[k+1], "fswatch"))   darktable.unmuted |= DT_DEBUG_FSWATCH; // fswatch module
        else if(!strcmp(argv[k+1], "camctl"))    darktable.unmuted |= DT_DEBUG_CAMCTL; // camera control module
        else if(!strcmp(argv[k+1], "perf"))      darktable.unmuted |= DT_DEBUG_PERF; // performance measurements
        else if(!strcmp(argv[k+1], "pwstorage")) darktable.unmuted |= DT_DEBUG_PWSTORAGE; // pwstorage module
        else if(!strcmp(argv[k+1], "opencl"))    darktable.unmuted |= DT_DEBUG_OPENCL;    // gpu accel via opencl
        else if(!strcmp(argv[k+1], "sql"))       darktable.unmuted |= DT_DEBUG_SQL; // SQLite3 queries
        else return usage(argv[0]);
        k ++;
      }
      else if(argv[k][1] == 't' && argc > k+1)
      {
        darktable.num_openmp_threads = CLAMP(atol(argv[k+1]), 1, 100);
        printf("[dt_init] using %d threads for openmp parallel sections\n", darktable.num_openmp_threads);
        k ++;
      }
    }
    else
    {
      image_to_load = argv[k];
    }
  }

#ifdef _OPENMP
  omp_set_num_threads(darktable.num_openmp_threads);
#endif

  g_type_init();

  /* check cput caps */
  dt_check_cpu(argc,argv);

  
#ifdef HAVE_GEGL
  (void)setenv("GEGL_PATH", DARKTABLE_DATADIR"/gegl:/usr/lib/gegl-0.0", 1);
  gegl_init(&argc, &argv);
#endif
  // thread-safe init:
  dt_exif_init();
  char datadir[1024];
  dt_get_user_config_dir (datadir,1024);
  char filename[1024];
  snprintf(filename, 1024, "%s/darktablerc", datadir);

  // Initialize the filesystem watcher
  darktable.fswatch=dt_fswatch_new();

#ifdef HAVE_GPHOTO2
  // Initialize the camera control
  darktable.camctl=dt_camctl_new();
#endif
  // has to go first for settings needed by all the others.
  darktable.conf = (dt_conf_t *)malloc(sizeof(dt_conf_t));
  dt_conf_init(darktable.conf, filename);

  // get max lighttable thumbnail size:
  darktable.thumbnail_size = CLAMPS(dt_conf_get_int("plugins/lighttable/thumbnail_size"), 160, 1300);
  // and make sure it can be mip-mapped all the way from mip4 to mip0
  darktable.thumbnail_size /= 16;
  darktable.thumbnail_size *= 16;

  // Initialize the password storage engine
  darktable.pwstorage=dt_pwstorage_new();

  // check and migrate database into new XDG structure
  char dbfilename[2048]= {0};
  gchar *conf_db = dt_conf_get_string("database");
  if (conf_db && conf_db[0] != '/')
  {
    const char *homedir = dt_get_home_dir();
    snprintf (dbfilename,2048,"%s/%s", homedir, conf_db);
    if (g_file_test (dbfilename, G_FILE_TEST_EXISTS))
    {
      fprintf(stderr, "[init] moving database into new XDG directory structure\n");
      // move database into place
      char destdbname[2048]= {0};
      snprintf(destdbname,2048,"%s/%s",datadir,"library.db");
      if(!g_file_test (destdbname,G_FILE_TEST_EXISTS))
      {
        rename(dbfilename,destdbname);
        dt_conf_set_string("database","library.db");
      }
    }
    g_free(conf_db);
  }

  // check and migrate the cachedir
  char cachefilename[2048]= {0};
  char cachedir[2048]= {0};
  gchar *conf_cache = dt_conf_get_string("cachefile");
  if (conf_cache && conf_cache[0] != '/')
  {
    const char *homedir = dt_get_home_dir();
    snprintf (cachefilename,2048,"%s/%s",homedir, conf_cache);
    if (g_file_test (cachefilename,G_FILE_TEST_EXISTS))
    {
      fprintf(stderr, "[init] moving cache into new XDG directory structure\n");
      char destcachename[2048]= {0};
      snprintf(destcachename,2048,"%s/%s",cachedir,"mipmaps");
      if(!g_file_test (destcachename,G_FILE_TEST_EXISTS))
      {
        rename(cachefilename,destcachename);
        dt_conf_set_string("cachefile","mipmaps");
      }
    }
    g_free(conf_cache);
  }

  gchar * dbname = NULL;
  if ( dbfilenameFromCommand == NULL )
  {
    dbname = dt_conf_get_string ("database");
    if(!dbname)               snprintf(dbfilename, 1024, "%s/library.db", datadir);
    else if(dbname[0] != '/') snprintf(dbfilename, 1024, "%s/%s", datadir, dbname);
    else                      snprintf(dbfilename, 1024, "%s", dbname);
  }
  else
  {
    snprintf(dbfilename, 1024, "%s", dbfilenameFromCommand);
    dbname = g_file_get_basename (g_file_new_for_path(dbfilenameFromCommand));
  }

  int load_cached = 1;
  // if db file does not exist, also don't load the cache.
  if(!g_file_test(dbfilename, G_FILE_TEST_IS_REGULAR)) load_cached = 0;
  if(sqlite3_open(dbfilename, &(darktable.db)))
  {
    fprintf(stderr, "[init] could not open database ");
    if(dbname) fprintf(stderr, "`%s'!\n", dbname);
    else       fprintf(stderr, "\n");
#ifndef HAVE_GCONF
    fprintf(stderr, "[init] maybe your %s/darktablerc is corrupt?\n",datadir);
    dt_get_datadir(dbfilename, 512);
    fprintf(stderr, "[init] try `cp %s/darktablerc %s/darktablerc'\n", dbfilename,datadir);
#else
    fprintf(stderr, "[init] check your /apps/darktable/database gconf entry!\n");
#endif
    sqlite3_close(darktable.db);
    if (dbname != NULL) g_free(dbname);
    return 1;
  }
  if (dbname != NULL) g_free(dbname);

  dt_pthread_mutex_init(&(darktable.db_insert), NULL);
  dt_pthread_mutex_init(&(darktable.plugin_threadsafe), NULL);

  darktable.control = (dt_control_t *)malloc(sizeof(dt_control_t));
  if(init_gui)
  {
    dt_control_init(darktable.control);
  }
  else
  {
    // this is in memory, so schema can't exist yet.
    if(!strcmp(dbfilename, ":memory:"))
    {
      dt_control_create_database_schema();
      dt_gui_presets_init(); // also init preset db schema.
    }
    darktable.control->running = 0;
    darktable.control->accels_global = NULL;
    darktable.control->accels_lighttable = NULL;
    darktable.control->accels_darkroom = NULL;
    darktable.control->accels_filmstrip = NULL;
    darktable.control->accels_capture = NULL;
    dt_pthread_mutex_init(&darktable.control->run_mutex, NULL);
  }

  // initialize collection query
  darktable.collection_listeners = NULL;
  darktable.collection = dt_collection_new(NULL);

  darktable.opencl = (dt_opencl_t *)malloc(sizeof(dt_opencl_t));
  dt_opencl_init(darktable.opencl, argc, argv);

  darktable.blendop = (dt_blendop_t *)malloc(sizeof(dt_blendop_t));
  dt_develop_blend_init(darktable.blendop);

  darktable.points = (dt_points_t *)malloc(sizeof(dt_points_t));
  dt_points_init(darktable.points, dt_get_num_threads());

  int thumbnails = dt_conf_get_int ("mipmap_cache_thumbnails");
  thumbnails = MIN(1000000, MAX(20, thumbnails));

  darktable.mipmap_cache = (dt_mipmap_cache_t *)malloc(sizeof(dt_mipmap_cache_t));
  dt_mipmap_cache_init(darktable.mipmap_cache, thumbnails);

  darktable.image_cache = (dt_image_cache_t *)malloc(sizeof(dt_image_cache_t));
  dt_image_cache_init(darktable.image_cache, MIN(10000, MAX(500, thumbnails)), load_cached);

  // The GUI must be initialized before the views, because the init()
  // functions of the views depend on darktable.control->accels_* to register
  // their keyboard accelerators

  if(init_gui)
  {
    darktable.gui = (dt_gui_gtk_t *)malloc(sizeof(dt_gui_gtk_t));
    if(dt_gui_gtk_init(darktable.gui, argc, argv)) return 1;
  }

  darktable.view_manager = (dt_view_manager_t *)malloc(sizeof(dt_view_manager_t));
  dt_view_manager_init(darktable.view_manager);

  // load the darkroom mode plugins once:
  dt_iop_load_modules_so();

  if(init_gui)
  {
    darktable.lib = (dt_lib_t *)malloc(sizeof(dt_lib_t));
    dt_lib_init(darktable.lib);

    dt_control_load_config(darktable.control);
    g_strlcpy(darktable.control->global_settings.dbname, filename, 512); // overwrite if relocated.

    darktable.imageio = (dt_imageio_t *)malloc(sizeof(dt_imageio_t));
    dt_imageio_init(darktable.imageio);
  }

  if(init_gui)
  {
    // Loading the keybindings
    char keyfile[1024];

    // First dump the default keymapping
    snprintf(keyfile, 1024, "%s/keyboardrc_default", datadir);
    gtk_accel_map_save(keyfile);

    // Removing extraneous semi-colons from the default keymap
    strip_semicolons_from_keymap(keyfile);

    // Then load any modified keys if available
    snprintf(keyfile, 1024, "%s/keyboardrc", datadir);
    if(g_file_test(keyfile, G_FILE_TEST_EXISTS))
      gtk_accel_map_load(keyfile);
    else
      gtk_accel_map_save(keyfile); // Save the default keymap if none is present
  }

  int id = 0;
  if(init_gui && image_to_load)
  {
    char* filename;
    if(g_str_has_prefix(image_to_load, "file://"))
      image_to_load += strlen("file://");
    if(g_path_is_absolute(image_to_load) == FALSE)
    {
      char* current_dir = g_get_current_dir();
      char* tmp_filename = g_build_filename(current_dir, image_to_load, NULL);
      filename = g_realpath(tmp_filename);
      if(filename == NULL)
      {
        dt_control_log(_("found strange path `%s'"), tmp_filename);
        g_free(current_dir);
        g_free(tmp_filename);
        g_free(filename);
        return 0;
      }
      g_free(current_dir);
      g_free(tmp_filename);
    }
    else
    {
      filename = g_strdup(image_to_load);
    }

    if(g_file_test(filename, G_FILE_TEST_IS_DIR))
    {
      // import a directory into a film roll
      unsigned int last_char = strlen(filename)-1;
      if(filename[last_char] == '/')
        filename[last_char] = '\0';
      id = dt_film_import(filename);
      if(id)
      {
        dt_film_open(id);
        dt_ctl_switch_mode_to(DT_LIBRARY);
      }
      else
      {
        dt_control_log(_("error loading directory `%s'"), filename);
      }
    }
    else
    {
      // import a single image
      gchar *directory = g_path_get_dirname((const gchar *)filename);
      dt_film_t film;
      const int filmid = dt_film_new(&film, directory);
      id = dt_image_import(filmid, filename, TRUE);
      g_free (directory);
      if(id)
      {
        dt_film_open(filmid);
        // make sure buffers are loaded (load full for testing)
        dt_image_t *img = dt_image_cache_get(id, 'r');
        dt_image_buffer_t buf = dt_image_get_blocking(img, DT_IMAGE_FULL, 'r');
        if(!buf)
        {
          id = 0;
          dt_image_cache_release(img, 'r');
          dt_control_log(_("file `%s' has unknown format!"), filename);
        }
        else
        {
          dt_image_release(img, DT_IMAGE_FULL, 'r');
          dt_image_cache_release(img, 'r');
          DT_CTL_SET_GLOBAL(lib_image_mouse_over_id, id);
          dt_ctl_switch_mode_to(DT_DEVELOP);
        }
      }
      else
      {
        dt_control_log(_("error loading file `%s'"), filename);
      }
    }
    g_free(filename);
  }
  if(init_gui && !id)
  {
    dt_ctl_switch_mode_to(DT_LIBRARY);
  }

  return 0;
}

void dt_cleanup()
{
  dt_ctl_switch_mode_to(DT_MODE_NONE);

  dt_control_write_config(darktable.control);
  dt_control_shutdown(darktable.control);

  dt_lib_cleanup(darktable.lib);
  free(darktable.lib);
  dt_view_manager_cleanup(darktable.view_manager);
  free(darktable.view_manager);
  dt_imageio_cleanup(darktable.imageio);
  free(darktable.imageio);
  dt_gui_gtk_cleanup(darktable.gui);
  free(darktable.gui);
  dt_image_cache_cleanup(darktable.image_cache);
  free(darktable.image_cache);
  dt_mipmap_cache_cleanup(darktable.mipmap_cache);
  free(darktable.mipmap_cache);
  dt_control_cleanup(darktable.control);
  free(darktable.control);
  dt_conf_cleanup(darktable.conf);
  free(darktable.conf);
  dt_points_cleanup(darktable.points);
  free(darktable.points);
  dt_iop_unload_modules_so();
  dt_opencl_cleanup(darktable.opencl);
  free(darktable.opencl);
#ifdef HAVE_GPHOTO2
  dt_camctl_destroy(darktable.camctl);
#endif
  dt_pwstorage_destroy(darktable.pwstorage);
  dt_fswatch_destroy(darktable.fswatch);

  sqlite3_close(darktable.db);
  dt_pthread_mutex_destroy(&(darktable.db_insert));
  dt_pthread_mutex_destroy(&(darktable.plugin_threadsafe));

  dt_exif_cleanup();
#ifdef HAVE_GEGL
  gegl_exit();
#endif
}

void dt_print(dt_debug_thread_t thread, const char *msg, ...)
{
  if(darktable.unmuted & thread)
  {
    va_list ap;
    va_start(ap, msg);
    vprintf(msg, ap);
    va_end(ap);
    fflush(stdout);
  }
}

void dt_gettime_t(char *datetime, time_t t)
{
  struct tm tt;
  (void)localtime_r(&t, &tt);
  strftime(datetime, 20, "%Y:%m:%d %H:%M:%S", &tt);
}

void dt_gettime(char *datetime)
{
  dt_gettime_t(datetime, time(NULL));
}

void *dt_alloc_align(size_t alignment, size_t size)
{
#if defined(__MACH__) || defined(__APPLE__) || (defined(__FreeBSD_version) && __FreeBSD_version < 700013)
  return malloc(size);
#elif defined(__WIN32__)
  return _aligned_malloc(size, alignment);
#else
  void *ptr = NULL;
  if(posix_memalign(&ptr, alignment, size)) return NULL;
  return ptr;
#endif
}

void
dt_get_user_config_dir (char *data, size_t bufsize)
{
  g_snprintf (data,bufsize,"%s/.config/darktable",dt_get_home_dir());
  if (g_file_test (data,G_FILE_TEST_EXISTS)==FALSE)
    g_mkdir_with_parents (data,0700);
}

void
dt_get_user_cache_dir (char *data, size_t bufsize)
{
  g_snprintf (data,bufsize,"%s/.cache/darktable",dt_get_home_dir());
  if (g_file_test (data,G_FILE_TEST_EXISTS)==FALSE)
    g_mkdir_with_parents (data,0700);
}


void
dt_get_user_local_dir (char *data, size_t bufsize)
{
  g_snprintf(data,bufsize,"%s/.local",dt_get_home_dir());
  if (g_file_test (data,G_FILE_TEST_EXISTS)==FALSE)
    g_mkdir_with_parents (data,0700);

}

void dt_get_plugindir(char *datadir, size_t bufsize)
{
#if defined(__MACH__) || defined(__APPLE__)
  gchar *curr = g_get_current_dir();
  int contains = 0;
  for(int k=0; darktable.progname[k] != 0; k++) if(darktable.progname[k] == '/')
    {
      contains = 1;
      break;
    }
  if(darktable.progname[0] == '/') // absolute path
    snprintf(datadir, bufsize, "%s", darktable.progname);
  else if(contains) // relative path
    snprintf(datadir, bufsize, "%s/%s", curr, darktable.progname);
  else
  {
    // no idea where we have been called. use compiled in path
    g_free(curr);
    snprintf(datadir, bufsize, "%s/darktable", DARKTABLE_LIBDIR);
    return;
  }
  size_t len = MIN(strlen(datadir), bufsize);
  char *t = datadir + len; // strip off bin/darktable
  for(; t>datadir && *t!='/'; t--);
  t--;
  if(*t == '.' && *(t-1) != '.')
  {
    for(; t>datadir && *t!='/'; t--);
    t--;
  }
  for(; t>datadir && *t!='/'; t--);
  g_strlcpy(t, "/lib/darktable", bufsize-(t-datadir));
  g_free(curr);
#else
  snprintf(datadir, bufsize, "%s/darktable", DARKTABLE_LIBDIR);
#endif
}

void dt_get_datadir(char *datadir, size_t bufsize)
{
#if defined(__MACH__) || defined(__APPLE__)
  gchar *curr = g_get_current_dir();
  int contains = 0;
  for(int k=0; darktable.progname[k] != 0; k++) if(darktable.progname[k] == '/')
    {
      contains = 1;
      break;
    }
  if(darktable.progname[0] == '/') // absolute path
    snprintf(datadir, bufsize, "%s", darktable.progname);
  else if(contains) // relative path
    snprintf(datadir, bufsize, "%s/%s", curr, darktable.progname);
  else
  {
    // no idea where we have been called. use compiled in path
    g_free(curr);
    snprintf(datadir, bufsize, "%s", DARKTABLE_DATADIR);
    return;
  }
  size_t len = MIN(strlen(datadir), bufsize);
  char *t = datadir + len; // strip off bin/darktable
  for(; t>datadir && *t!='/'; t--);
  t--;
  if(*t == '.' && *(t-1) != '.')
  {
    for(; t>datadir && *t!='/'; t--);
    t--;
  }
  for(; t>datadir && *t!='/'; t--);
  g_strlcpy(t, "/share/darktable", bufsize-(t-datadir));
  g_free(curr);
#else
  snprintf(datadir, bufsize, "%s", DARKTABLE_DATADIR);
#endif
}

void dt_show_times(const dt_times_t *start, const char *prefix, const char *suffix, ...)
{
  dt_times_t end;
  char buf[120];		/* Arbitrary size, should be lots big enough for everything used in DT */
  int i;

  /* Skip all the calculations an everything if -d perf isn't on */
  if (darktable.unmuted & DT_DEBUG_PERF)
  {
    dt_get_times(&end);
    i = sprintf(buf, "%s took %.3f secs (%.3f CPU)", prefix, end.clock - start->clock, end.user - start->user);
    if (suffix != NULL)
    {
      va_list ap;
      va_start(ap, suffix);
      buf[i++] = ' ';
      vsnprintf(buf + i, sizeof buf - i, suffix, ap);
      va_end(ap);
    }
    dt_print(DT_DEBUG_PERF, "%s\n", buf);
  }
}

// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
