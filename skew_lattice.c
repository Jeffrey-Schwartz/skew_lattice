/*
 *  @(#) $Id: skew_lattice.c 2014-05-08 $
 *  Copyright (C) 2014 Jeffrey J. Schwartz.
 *  E-mail: schwartz@physics.ucla.edu
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
 */

/*
 *  This program skews SPM images in an attempt to compensate for
 *  lateral drift during imaging. Skewing the images by sequential
 *  lateral translations of subsequent rows/cols in the image is 
 *  used to regularize the lattice to its known parameters. Angles
 *  between lattice features may be measured in the program to aid
 *  the user to determining the optimal skew amount.
 */

#include "config.h"
#include <string.h>
#include <gtk/gtk.h>
#include <app/gwyapp.h>
#include <app/gwymoduleutils.h>
#include <libprocess/stats.h>
#include <libprocess/filters.h>
#include <libprocess/inttrans.h>
#include <libprocess/datafield.h>
#include <libprocess/gwyprocess.h>
#include <libgwyddion/gwymath.h>
#include <libgwydgets/gwydataview.h>
#include <libgwydgets/gwydgetutils.h>
#include <libgwydgets/gwynullstore.h>
#include <libgwydgets/gwylayer-basic.h>
#include <libgwydgets/gwyradiobuttons.h>
#include <libgwymodule/gwymodule-process.h>

#define skew_lattice_RUN_MODES (GWY_RUN_INTERACTIVE)
#define PI 3.14159265358979323846

typedef struct _GwyToolLevel3      GwyToolLevel3;

typedef struct _GwyToolLevel3Class GwyToolLevel3Class;

struct _GwyToolLevel3
{
    GwyPlainTool parent_instance;
    GtkTreeView *treeview;
    GtkTreeModel *model;
    GtkObject *radius;
    gint32 rpx;
    GtkWidget *instant_apply;
    GtkWidget *set_zero;
    GtkWidget *apply;
    GType layer_type_point;
};

enum
{
    COLUMN_I,
    COLUMN_X,
    COLUMN_Y,
    COLUMN_Z,
    NCOLUMNS
};

enum
{
    PREVIEW_SIZE = 512
};

typedef enum {
    IMAGE_DATA,
    IMAGE_FFT,
    IMAGE_CORRECTED,
    IMAGE_FFT_CORRECTED,
} ImageMode;

typedef enum {
    ZOOM_1 = 1,
    ZOOM_2 = 2,
} ZoomMode;

typedef enum {
    HORIZONTAL,
    VERTICAL,
} ShiftMode;

typedef struct {
    gdouble lower;
    gdouble upper;
    gfloat Xskew;
    gfloat Yskew;
    gdouble angle1;
    gdouble angle2;
    ImageMode image_mode;
    ZoomMode zoom_mode;
    gint copy_row_start;
    gint copy_col_start;
    gdouble background_fill;
    gboolean background;
    gint newxres;
    gint newyres;
} ThresholdArgs;

typedef struct {
    gdouble min, max;
} ThresholdRanges;

typedef struct {
    ThresholdArgs *args;
    ThresholdRanges *ranges;
    GtkWidget *dialog;
    GtkWidget *view;
    GtkWidget *lower;
    GtkWidget *upper;
    GtkWidget *hskewtxt;
    GtkWidget *vskewtxt;
    GwyContainer *mydata;
    GwyContainer *container;
    GwyDataField *dfield;
    GwyDataField *image;
    GwyDataField *corr_image;
    GwyDataField *corr_fft;
    GwyDataField *disp_data;
    gint id;
    GwySelection *selection;
    GwySIValueFormat *original_XY_Format;
    GwySIValueFormat *XY_Format;
    GwySIValueFormat *Z_Format;
    GwySIUnit *Image_XY_Units;
    GwySIUnit *Image_Z_Units;
    GwyToolLevel3 *tool;
    GSList *image_mode_radios;
    GSList *zoom_mode_radios;
    GtkObject *skew_Xadjust;
    GtkWidget *skew_Xslider;
    GtkObject *skew_Yadjust;
    GtkWidget *skew_Yslider;
    GtkWidget *Angle1;
    GtkWidget *Angle2;
    gdouble p[4][3];
    GwyVectorLayer *vlayer;
} ThresholdControls;

static gboolean module_register             (void);

static void     skew_lattice                 (GwyContainer *data, GwyRunType run);
static void     perform_fft                 (GwyDataField *dfield,
                                                GwyContainer *data);
static void     selection_changed           (ThresholdControls *controls);
static void     clear_points                (ThresholdControls *controls);
static void     peak_find                   (ThresholdControls *controls,
                                                gdouble *point, guint idx);
static void     skew_do                     (ThresholdControls *controls);
static void     skew_create_output          (GwyContainer *data, 
                                                GwyDataField *dfield,
                                                ThresholdControls *controls);
static void     skew_lattice_dialog             (ThresholdControls *controls,
                                            ThresholdRanges *ranges,
                                            GwyContainer *data,
                                            GwyDataField *dfield,
                                            gint id);
static void     threshold_set_to_full_range(ThresholdControls *controls);
static void     threshold_lower_changed    (ThresholdControls *controls);
static void     threshold_upper_changed    (ThresholdControls *controls);
static void     preview                    (ThresholdControls *controls);
static void     threshold_do               (ThresholdArgs *args,
                                            GwyDataField *dfield);
static void     threshold_load_args        (ThresholdControls *controls);
static void     threshold_save_args        (ThresholdControls *controls);
static void     zoom_mode_changed          (GtkToggleButton *button,
                                            ThresholdControls *controls);
static void     image_mode_changed         (GtkToggleButton *button,
                                            ThresholdControls *controls);
static void     gwy_tool_level3_render_cell    (GtkCellLayout *layout,
                            GtkCellRenderer *renderer, GtkTreeModel *model,
                            GtkTreeIter *iter, gpointer user_data);
static void     gwy_tool_level3_radius_changed(GwyToolLevel3 *tool);
static void     fft_postprocess            (GwyDataField *dfield);
static void     set_dfield_modulus         (GwyDataField *re, GwyDataField *im,
                                                GwyDataField *target);
static void radio_buttons_attach_to_table  (GSList *group,
                                                GtkTable *table, gint row);
static void     skew_update_angles      (ThresholdControls *controls);
static void     get_angles              (ThresholdControls *controls);
static void     reFind_Peaks            (ThresholdControls *controls);
static void     zoom_adjust_peaks       (ThresholdControls *controls);
static void     skew_Xadjusted          (ThresholdControls *controls);
static void     skew_Yadjusted          (ThresholdControls *controls);
static void     skew_process            (ThresholdControls *controls);
static void     reset_Xskew             (ThresholdControls *controls);
static void     reset_Yskew             (ThresholdControls *controls);
static void     hskew_changed           (ThresholdControls *controls);
static void     vskew_changed           (ThresholdControls *controls);
static gdouble  matrix_det              (const gdouble *m);
static void     invert_matrix           (gdouble *dest, const gdouble *src);
static gdouble  deg2rad                 (const gdouble deg);
static void     mult_3matrix            (gdouble *dest, const gdouble *mat1,
                                        const gdouble *mat2);
static void     affine                  (GwyDataField *source,
                                        GwyDataField *dest,
                                        const gdouble *invtrans,
                                        GwyInterpolationType interp,
                                        gdouble fill_value);

static const ThresholdArgs threshold_defaults = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3, 1, 0, 0, 0.0, FALSE, 0, 0
};


static GwyModuleInfo module_info = {
    GWY_MODULE_ABI_VERSION, &module_register,
    N_("Tool to correct for lateral drift during scanning probe imaging; "
    "skews image to obtain regular lattice shape."),
    "Jeffrey J. Schwartz <schwartz@physics.ucla.edu>",
    "1.0",
    "Jeffrey J. Schwartz",
    "May 2014",
};

GWY_MODULE_QUERY(module_info)

static gboolean
module_register(void)
{
    gwy_process_func_register("skew_lattice",
                (GwyProcessFunc)&skew_lattice,
                N_("/_Correct Data/_Skew Lattice"),
                NULL, skew_lattice_RUN_MODES, GWY_MENU_FLAG_DATA,
                N_("Skews image to form regular lattice"));
    return TRUE;
}

static void
skew_lattice(GwyContainer *data, GwyRunType run)
{
    ThresholdControls controls;
    ThresholdArgs args;
    ThresholdRanges ranges;
    GwyDataField *dfield;
    GQuark quark;
    gint id;
    GwyToolLevel3 tool;
    tool.rpx = 3;
    args = threshold_defaults;
    controls.args = &args;
    controls.tool = &tool;
    g_return_if_fail(run & skew_lattice_RUN_MODES);
    gwy_app_data_browser_get_current(GWY_APP_DATA_FIELD, &dfield,
                                     GWY_APP_DATA_FIELD_ID, &id,
                                     GWY_APP_DATA_FIELD_KEY, &quark, 0);
    g_return_if_fail(dfield);
    if (run == GWY_RUN_INTERACTIVE)
    {
        skew_lattice_dialog(&controls, &ranges, data,
            gwy_data_field_duplicate(dfield), id);
        gwy_data_field_data_changed(dfield);
    }
}

static void
threshold_format_value(ThresholdControls *controls,
                       GtkEntry *entry, gdouble value)
{
    gchar *s;
    s = g_strdup_printf("%.*f",
                        controls->original_XY_Format->precision+1,
                        value/controls->original_XY_Format->magnitude);
    gtk_entry_set_text(GTK_ENTRY(entry), s);
    g_free(s);
}

static GtkWidget*
threshold_entry_attach(ThresholdControls *controls,
                       GtkTable *table, gint row,
                       gdouble value, const gchar *name)
{
    GtkWidget *label, *entry;
    label = gtk_label_new_with_mnemonic(name);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    entry = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(entry, TRUE);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 8);
    threshold_format_value(controls, GTK_ENTRY(entry), value);
    gtk_table_attach(table, entry, 1, 3, row, row+1, GTK_FILL, 0, 0, 0);
    label = gtk_label_new(controls->original_XY_Format->units);
    gtk_label_set_markup(GTK_LABEL(label),
                                    controls->original_XY_Format->units);
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 3, 4, row, row+1, GTK_FILL, 0, 0, 0);
    return entry;
}

static void
skew_lattice_dialog(ThresholdControls *controls, ThresholdRanges *ranges,
                 GwyContainer *data, GwyDataField *dfield, gint id)
{
    GtkWidget *dialog, *hbox, *button, *label;
    GtkTable *table;
    GwyVectorLayer *vlayer;
    gint response, row;
    GwyPixmapLayer *layer;
    controls->image = gwy_data_field_duplicate(dfield);
    controls->corr_image = gwy_data_field_duplicate(controls->image);
    controls->container = data;
    controls->id = id;    
    controls->ranges = ranges;
    controls->dfield = dfield;
    controls->disp_data = gwy_data_field_new_alike(dfield, TRUE);
    controls->original_XY_Format = gwy_data_field_get_value_format_xy
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    controls->Image_XY_Units = gwy_data_field_get_si_unit_xy(controls->image);
    controls->Image_Z_Units = gwy_data_field_get_si_unit_z(controls->image);
    controls->mydata = gwy_container_new();
    perform_fft(controls->dfield, controls->mydata);
    controls->corr_fft = gwy_data_field_duplicate(controls->dfield);
    gwy_data_field_get_min_max(dfield, &ranges->min, &ranges->max);
    controls->XY_Format = gwy_data_field_get_value_format_xy
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    controls->Z_Format = gwy_data_field_get_value_format_z
                            (dfield, GWY_SI_UNIT_FORMAT_MARKUP, NULL);
    dialog = gtk_dialog_new_with_buttons(_("Skew Lattice"), NULL, 0,
                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK, GTK_RESPONSE_OK, NULL);
    gtk_dialog_set_default_response(GTK_DIALOG(dialog), GTK_RESPONSE_OK);
    controls->dialog = dialog;
    hbox = gtk_hbox_new(FALSE, 2);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox,
                       FALSE, FALSE, 4);
    table = GTK_TABLE(gtk_table_new(4, 4, FALSE));
    gtk_table_set_row_spacings(table, 2);
    gtk_table_set_col_spacings(table, 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);
    gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(table), TRUE, TRUE, 4);
    label = gtk_label_new("Data Display");
    gtk_label_set_markup(GTK_LABEL(label),
        "<b>Data Display</b>\n(FFT: Modulus, Hanning window, subtract mean)");
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment(GTK_MISC(label), 0.5, 0.5);
    gtk_table_attach(table, label, 0, 4, 0, 1, GTK_FILL, 0, 0, 0);
    gwy_app_sync_data_items(data, controls->mydata, id, 0, FALSE,
                GWY_DATA_ITEM_PALETTE, GWY_DATA_ITEM_MASK_COLOR,
                GWY_DATA_ITEM_RANGE, GWY_DATA_ITEM_REAL_SQUARE, 0);
    gwy_container_set_object_by_name(controls->mydata, "/0/data", dfield);
    controls->view = gwy_data_view_new(controls->mydata);
    layer = gwy_layer_basic_new();
    g_object_set(layer, "data-key", "/0/data",
                 "gradient-key", "/0/base/palette",
                 "range-type-key", "/0/base/range-type",
                 "min-max-key", "/0/base", NULL);
    gwy_data_view_set_data_prefix(GWY_DATA_VIEW(controls->view), "/0/data");
    gwy_data_view_set_base_layer(GWY_DATA_VIEW(controls->view), layer);
    gwy_set_data_preview_size(GWY_DATA_VIEW(controls->view), PREVIEW_SIZE);
    vlayer = g_object_new(g_type_from_name("GwyLayerPoint"),
                  "selection-key", "/0/select/point", NULL);
    controls->vlayer = vlayer;
    gwy_data_view_set_top_layer(GWY_DATA_VIEW(controls->view), vlayer);
    controls->selection = gwy_vector_layer_ensure_selection(vlayer);
    gwy_selection_set_max_objects(controls->selection, 4);
    g_signal_connect_swapped(controls->selection, "changed",
                         G_CALLBACK(selection_changed), controls);
    gtk_table_attach(table, controls->view, 0, 4, 1, 2, GTK_FILL, 0, 0, 0);
    label = gtk_label_new("Select four sequential peaks "
                          "in the first ring around center");
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    gtk_misc_set_alignment(GTK_MISC(label), 0.5, 0.5);
    gtk_table_attach(table, label, 0, 4, 2, 3, GTK_FILL, 0, 0, 0);
    controls->Angle1 = gtk_label_new("Angle 123:");
    gtk_label_set_markup(GTK_LABEL(controls->Angle1), "<b>Angle 123:</b>");
    gtk_label_set_width_chars (GTK_LABEL(controls->Angle1), 15);
    gtk_misc_set_alignment(GTK_MISC(controls->Angle1), 0.0, 0.0);
    gtk_table_attach(table, controls->Angle1, 2, 3,
                                            3, 4, GTK_FILL, 0, 0, 0);
    controls->Angle2 = gtk_label_new("Angle 234:");
    gtk_label_set_markup(GTK_LABEL(controls->Angle2), "<b>Angle 234:</b>");
    gtk_label_set_width_chars (GTK_LABEL(controls->Angle2), 15);
    gtk_misc_set_alignment(GTK_MISC(controls->Angle2), 0.0, 1.0);
    gtk_table_attach(table, controls->Angle2, 3, 4,
                                            3, 4, GTK_FILL, 0, 0, 0);
    table = GTK_TABLE(gtk_table_new(7, 4, FALSE));
    gtk_table_set_row_spacings(table, 2);
    gtk_table_set_col_spacings(table, 6);
    gtk_container_set_border_width(GTK_CONTAINER(table), 4);
    gtk_box_pack_start(GTK_BOX(hbox), GTK_WIDGET(table), TRUE, TRUE, 4);
    row = 0;
    label = gtk_label_new("Display Zoom: ");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Zoom:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
    gtk_table_attach(table, label, 0, 1, row, row+1, GTK_FILL, 0, 0, 0);
    row++;    
    controls->zoom_mode_radios
        = gwy_radio_buttons_createl(G_CALLBACK(zoom_mode_changed), controls,
                                    controls->args->zoom_mode,
                                    _("×1"), ZOOM_1,
                                    _("×2"), ZOOM_2,
                                    NULL);
    radio_buttons_attach_to_table(controls->zoom_mode_radios, table, row);
    row++;
    label = gtk_label_new("Specify intensity range:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Specify intensity range:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 7, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    controls->lower = threshold_entry_attach(controls, table,
                             row, controls->args->lower, _("_Lower:"));
    g_signal_connect_swapped(controls->lower, "activate",
                             G_CALLBACK(threshold_lower_changed), controls);
    row++;
    controls->upper = threshold_entry_attach(controls, table, row,
                            controls->args->upper, _("_Upper:"));
    g_signal_connect_swapped(controls->upper, "activate",
                            G_CALLBACK(threshold_upper_changed), controls);
    row++;
    button = gtk_button_new_with_mnemonic(_("Set to _Full Range"));
    gtk_table_attach(table, button, 0, 4, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(threshold_set_to_full_range),
                             controls);
    row++;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 20);
    label = gtk_label_new("Peak Positions:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Peak positions:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 4, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    GtkTreeViewColumn *column;
    GtkCellRenderer *renderer;
    GwyNullStore *store;
    store = gwy_null_store_new(4);
    controls->tool->model = GTK_TREE_MODEL(store);
    controls->tool->treeview =
        GTK_TREE_VIEW(gtk_tree_view_new_with_model(controls->tool->model));
    gchar *XUnits, *YUnits, *ZUnits;
    XUnits = g_strdup_printf("<b>x</b> [%s]", controls->XY_Format->units);
    YUnits = g_strdup_printf("<b>y</b> [%s]", controls->XY_Format->units);
    ZUnits = g_strdup_printf("<b>value</b> [%s]", controls->Z_Format->units);
    guint i;
    for (i = 0; i < NCOLUMNS; i++) {
        column = gtk_tree_view_column_new();
        g_object_set_data(G_OBJECT(column), "id", GUINT_TO_POINTER(i));
        renderer = gtk_cell_renderer_text_new();
        g_object_set(renderer, "xalign", 1.0, NULL);
        gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(column), renderer, TRUE);
        gtk_cell_layout_set_cell_data_func(GTK_CELL_LAYOUT(column), renderer,
                                            gwy_tool_level3_render_cell,
                                            controls, NULL);
        label = gtk_label_new(NULL);
        switch (i)
        {
            case 0:
                gtk_label_set_markup(GTK_LABEL(label), "<b>n</b>");
                break;
            case 1:
                gtk_label_set_markup(GTK_LABEL(label), XUnits);
                break;
            case 2:
                gtk_label_set_markup(GTK_LABEL(label), YUnits);
                break;
            case 3:
                gtk_label_set_markup(GTK_LABEL(label), ZUnits);
                break;
        }
        gtk_tree_view_column_set_widget(column, label);
        gtk_widget_show(label);
        gtk_tree_view_append_column(controls->tool->treeview, column);
    }
    gtk_table_attach(table, GTK_WIDGET(controls->tool->treeview),
            0, 4, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    g_free(XUnits);
    g_free(YUnits);
    g_free(ZUnits);
    button = gtk_button_new_with_mnemonic(_("Clear Points"));
    gtk_table_attach(table, button, 0, 4, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(clear_points), controls);
    row++;
    controls->tool->radius = gtk_adjustment_new(controls->tool->rpx,
                                                            0, 10, 1, 5, 0);
    gwy_table_attach_spinbutton(GTK_WIDGET(table), row, 
                _("Peak search radius:"),
                "px", controls->tool->radius);
    g_signal_connect_swapped(controls->tool->radius, "value-changed",
                 G_CALLBACK(gwy_tool_level3_radius_changed), controls->tool);
    row++;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 10);
    label = gtk_label_new("Display Mode:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Display Mode:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 5, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    controls->image_mode_radios
        = gwy_radio_buttons_createl(G_CALLBACK(image_mode_changed), controls,
                                    controls->args->image_mode,
                                    _("Image"), IMAGE_DATA,
                                    _("Image FFT"), IMAGE_FFT,
                                    _("Skewed Image"), IMAGE_CORRECTED,
                                    _("Skewed FFT"), IMAGE_FFT_CORRECTED,
                                    NULL);
    radio_buttons_attach_to_table(controls->image_mode_radios, table, row);
    row += 2;
    gtk_table_set_row_spacing(GTK_TABLE(table), row-1, 10);
    label = gtk_label_new("Horizontal Skew:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Horizontal Skew:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    button = gtk_button_new_with_mnemonic(_("Reset X Skew"));
    gtk_table_attach(table, button, 3, 5, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(reset_Xskew), controls);
    row++;
    controls->skew_Xadjust = gtk_adjustment_new(0, -30, 30, 1, 1, 0);
    controls->skew_Xslider = gtk_hscale_new(
                                (GtkAdjustment*)controls->skew_Xadjust);
    gtk_table_attach(table, controls->skew_Xslider, 0, 3,
                            row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(controls->skew_Xadjust, "value-changed",
                         G_CALLBACK(skew_Xadjusted), controls);
    controls->hskewtxt = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(controls->hskewtxt, TRUE);
    gtk_entry_set_width_chars(GTK_ENTRY(controls->hskewtxt), 5);
    gtk_entry_set_text(GTK_ENTRY(controls->hskewtxt), "0.0");
    gtk_table_attach(table, controls->hskewtxt, 3, 4,
                            row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(controls->hskewtxt, "activate",
                            G_CALLBACK(hskew_changed), controls);
    label = gtk_label_new("deg");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
    gtk_table_attach(table, label, 4, 5, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    label = gtk_label_new("Vertical Skew:");
    gtk_label_set_markup(GTK_LABEL(label), "<b>Vertical Skew:</b>");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach(table, label, 0, 3, row, row+1, GTK_FILL, 0, 0, 0);
    button = gtk_button_new_with_mnemonic(_("Reset Y Skew"));
    gtk_table_attach(table, button, 3, 5, row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(button, "clicked",
                             G_CALLBACK(reset_Yskew), controls);
    row++;
    controls->skew_Yadjust = gtk_adjustment_new(0, -30, 30, 1, 1, 0);
    controls->skew_Yslider = gtk_hscale_new(
                                (GtkAdjustment*)controls->skew_Yadjust);
    gtk_table_attach(table, controls->skew_Yslider, 0, 3,
                            row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(controls->skew_Yadjust, "value-changed",
                         G_CALLBACK(skew_Yadjusted), controls);
    controls->vskewtxt = gtk_entry_new();
    gwy_widget_set_activate_on_unfocus(controls->vskewtxt, TRUE);
    gtk_entry_set_width_chars(GTK_ENTRY(controls->vskewtxt), 5);
    gtk_entry_set_text(GTK_ENTRY(controls->vskewtxt), "0.0");
    gtk_table_attach(table, controls->vskewtxt, 3, 4,
                            row, row+1, GTK_FILL, 0, 0, 0);
    g_signal_connect_swapped(controls->vskewtxt, "activate",
                             G_CALLBACK(vskew_changed), controls);
    label = gtk_label_new("deg");
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.0);
    gtk_table_attach(table, label, 4, 5, row, row+1, GTK_FILL, 0, 0, 0);
    row++;
    threshold_load_args(controls);
    skew_process(controls);
    preview(controls);
    gtk_widget_show_all(dialog);
    do
    {
        response = gtk_dialog_run(GTK_DIALOG(dialog));
        switch (response)
        {
            case GTK_RESPONSE_CANCEL:
            case GTK_RESPONSE_DELETE_EVENT:
                gtk_widget_destroy(dialog);
            case GTK_RESPONSE_NONE:
                g_object_unref(controls->mydata);
                gwy_si_unit_value_format_free(controls->XY_Format);
                gwy_si_unit_value_format_free(controls->Z_Format);
                threshold_save_args(controls);
                return;
                break;
            case GTK_RESPONSE_OK:
                break;
            default:
                g_assert_not_reached();
                break;
        }
    } while (response != GTK_RESPONSE_OK);
    threshold_save_args(controls);
    skew_do(controls);
    gtk_widget_destroy(dialog);
    g_object_unref(controls->mydata);
    gwy_si_unit_value_format_free(controls->original_XY_Format);
    gwy_si_unit_value_format_free(controls->XY_Format);
    gwy_si_unit_value_format_free(controls->Z_Format);
}

static void
radio_buttons_attach_to_table(GSList *group, 
                GtkTable *table, gint row)
{
    g_return_val_if_fail(GTK_IS_TABLE(table), row);
    while (group)
    {
        gtk_table_attach(table, GTK_WIDGET(group->data),
                         0, 2, row, row + 1,
                         GTK_EXPAND | GTK_FILL, 0, 0, 0);
        group = g_slist_next(group);
        gtk_table_attach(table, GTK_WIDGET(group->data),
                         3, 5, row, row + 1,
                         GTK_EXPAND | GTK_FILL, 0, 0, 0);
        row++;
        group = g_slist_next(group);
    }
}

static void
threshold_set_to_range(ThresholdControls *controls,
                       gdouble lower, gdouble upper)
{
    threshold_format_value(controls, GTK_ENTRY(controls->lower), lower);
    gtk_widget_activate(controls->lower);
    threshold_format_value(controls, GTK_ENTRY(controls->upper), upper);
    gtk_widget_activate(controls->upper);
    preview(controls);
}

static void
threshold_set_to_full_range(ThresholdControls *controls)
{
    threshold_set_to_range(controls,
           controls->ranges->min, controls->ranges->max);
}

static void
threshold_lower_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->lower));
    gdouble num =
        g_strtod(value, NULL) * controls->original_XY_Format->magnitude;
    if (num >= controls->ranges->min && num <= controls->ranges->max)
        controls->args->lower = num;
    else
    {
        if (num < controls->ranges->min)
            controls->args->lower = controls->ranges->min;
        else if (num > controls->ranges->max)
            controls->args->lower = controls->ranges->max; 
    }
    threshold_format_value(controls,
        GTK_ENTRY(controls->lower), controls->args->lower);
    threshold_save_args(controls);
    preview(controls);
}

static void
threshold_upper_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->upper));
    gdouble num =
        g_strtod(value, NULL) * controls->original_XY_Format->magnitude;
    if (num >= controls->ranges->min && num <= controls->ranges->max)
        controls->args->upper = num;
    else
    {
        if (num < controls->ranges->min)
            controls->args->upper = controls->ranges->min;
        else if (num > controls->ranges->max)
            controls->args->upper = controls->ranges->max;
    }
    threshold_format_value(controls, GTK_ENTRY(controls->upper),
            controls->args->upper);
    threshold_save_args(controls);
    preview(controls);
}

static void
preview(ThresholdControls *controls)
{
    gint Xres, Yres;
    gdouble Xreal, Yreal;
    gdouble Xoff, Yoff;
    GwySIUnit *XY_Units;
    GwySIUnit *Z_Units;
    Xres = gwy_data_field_get_xres(controls->disp_data);
    Yres = gwy_data_field_get_yres(controls->disp_data);
    Xreal = gwy_data_field_get_xreal(controls->disp_data);
    Yreal = gwy_data_field_get_yreal(controls->disp_data);
    Xoff = gwy_data_field_get_xoffset(controls->disp_data);
    Yoff = gwy_data_field_get_yoffset(controls->disp_data);
    XY_Units = gwy_data_field_get_si_unit_xy(controls->disp_data);
    Z_Units = gwy_data_field_get_si_unit_z(controls->disp_data);
    switch (controls->args->image_mode)
    {
        case IMAGE_DATA:
            gwy_data_field_resample(controls->image,
                        Xres, Yres, GWY_INTERPOLATION_BILINEAR);
            gwy_data_field_copy(controls->image,
                                    controls->disp_data, TRUE);
            Xreal = gwy_data_field_get_xreal(controls->image);
            Yreal = gwy_data_field_get_yreal(controls->image);
            XY_Units = gwy_data_field_get_si_unit_xy(controls->image);
            Z_Units = gwy_data_field_get_si_unit_z(controls->image);
            Xoff = gwy_data_field_get_xoffset(controls->image);
            Yoff = gwy_data_field_get_yoffset(controls->image);
            break;
        case IMAGE_FFT:
            gwy_data_field_resample(controls->dfield,
                        Xres, Yres, GWY_INTERPOLATION_BILINEAR);
            gwy_data_field_copy(controls->dfield,
                                    controls->disp_data, TRUE);
            Xreal = gwy_data_field_get_xreal(controls->dfield);
            Yreal = gwy_data_field_get_yreal(controls->dfield);
            XY_Units = gwy_data_field_get_si_unit_xy(controls->dfield);
            Z_Units = gwy_data_field_get_si_unit_z(controls->dfield);
            Xoff = gwy_data_field_get_xoffset(controls->dfield);
            Yoff = gwy_data_field_get_yoffset(controls->dfield);
            break;
        case IMAGE_CORRECTED:
            gwy_data_field_resample(controls->corr_image,
                        Xres, Yres, GWY_INTERPOLATION_BILINEAR);
            gwy_data_field_copy(controls->corr_image,
                                    controls->disp_data, TRUE);
            Xreal = gwy_data_field_get_xreal(controls->corr_image);
            Yreal = gwy_data_field_get_yreal(controls->corr_image);
            XY_Units = gwy_data_field_get_si_unit_xy(controls->corr_image);
            Z_Units = gwy_data_field_get_si_unit_z(controls->corr_image);
            Xoff = gwy_data_field_get_xoffset(controls->corr_image);
            Yoff = gwy_data_field_get_yoffset(controls->corr_image);
            break;
        case IMAGE_FFT_CORRECTED:
            gwy_data_field_resample(controls->corr_fft,
                        Xres, Yres, GWY_INTERPOLATION_BILINEAR);
            gwy_data_field_copy(controls->corr_fft,
                                    controls->disp_data, TRUE);
            Xreal = gwy_data_field_get_xreal(controls->corr_fft);
            Yreal = gwy_data_field_get_yreal(controls->corr_fft);
            XY_Units = gwy_data_field_get_si_unit_xy(controls->corr_fft);
            Z_Units = gwy_data_field_get_si_unit_z(controls->corr_fft);
            Xoff = gwy_data_field_get_xoffset(controls->corr_fft);
            Yoff = gwy_data_field_get_yoffset(controls->corr_fft);
            break;
    }
    ZoomMode zoom = controls->args->zoom_mode;
    if (zoom != ZOOM_1)
    {
        guint width = (Xres/controls->args->zoom_mode) | 1;
        guint height = (Yres/controls->args->zoom_mode) | 1;
        GwyDataField *temp = gwy_data_field_area_extract(controls->disp_data,
                                          (Xres - width)/2, (Yres - height)/2,
                                          width, height);
        gwy_data_field_resample(temp, Xres, Yres, GWY_INTERPOLATION_BILINEAR);
        g_object_unref(controls->disp_data);
        controls->disp_data = gwy_data_field_duplicate(temp);
        g_object_unref(temp);
    }
    gwy_data_field_set_xreal(controls->disp_data, Xreal/zoom);
    gwy_data_field_set_yreal(controls->disp_data, Yreal/zoom);
    gwy_data_field_set_xoffset(controls->disp_data, Xoff/zoom);
    gwy_data_field_set_yoffset(controls->disp_data, Yoff/zoom);
    gwy_data_field_set_si_unit_xy(controls->disp_data, XY_Units);
    gwy_data_field_set_si_unit_z(controls->disp_data, Z_Units);
    gwy_container_set_object_by_name(controls->mydata,
                                    "/0/data", controls->disp_data);
    gwy_data_field_get_min_max(controls->disp_data, &controls->ranges->min,
                                    &controls->ranges->max);
    threshold_do(controls->args, controls->disp_data);
    gwy_set_data_preview_size(GWY_DATA_VIEW(controls->view), PREVIEW_SIZE);
}

static void
peak_find(ThresholdControls *controls, gdouble *point, guint idx)
{
    GwyDataField *dfield = controls->disp_data;
    gint i, j, low_i, high_i, low_j, high_j;
    gdouble temp_i, temp_j, temp_z;
    gint col = gwy_data_field_rtoj(dfield, point[0]);
    gint row = gwy_data_field_rtoi(dfield, point[1]);
    temp_i = col;
    temp_j = row;    
    temp_z = gwy_data_field_get_val(dfield, col, row);
    gint32 r = controls->tool->rpx;
    gint Xres, Yres;
    Xres = gwy_data_field_get_xres(dfield);
    Yres = gwy_data_field_get_yres(dfield);
    low_i = col - r;
    high_i = col + r;
    if (low_i < 0)
        low_i = 0;
    if (high_i > Xres)
        high_i = Xres;
    low_j = row - r;
    high_j = row + r;
    if (low_j < 0)
        low_j = 0;
    if (high_j > Yres)
        high_j = Yres;
    for (i = low_i; i < high_i; i++)
    {
        for (j = low_j; j < high_j; j++)
        {
            if (gwy_data_field_get_val(dfield, i, j) > temp_z)
            {
                temp_i = i;
                temp_j = j;
                temp_z = gwy_data_field_get_val(dfield, i, j);
            }
        }
    }
    controls->p[idx][0] = gwy_data_field_jtor(dfield, temp_i)
                + gwy_data_field_get_xoffset(dfield);
    controls->p[idx][1] = gwy_data_field_itor(dfield, temp_j)
                + gwy_data_field_get_yoffset(dfield);
    controls->p[idx][2] = temp_z;
    if ((row - temp_j) != 0 || (col - temp_i) != 0)
    {
        point[0] = gwy_data_field_jtor(dfield, temp_i);
        point[1] = gwy_data_field_itor(dfield, temp_j);
        gwy_selection_set_object(controls->selection, idx, point);
    }
}

static void
reFind_Peaks(ThresholdControls *controls)
{
    int i, num = 0;
    for (i = 0; i < 4; i++)
    {
        double point[2];
        if (gwy_selection_get_object(controls->selection, i, point))
        {
            gdouble xoff, yoff;
            peak_find(controls, point, i);
            xoff = gwy_data_field_get_xoffset(controls->disp_data);
            yoff = gwy_data_field_get_yoffset(controls->disp_data);
            point[0] = controls->p[num][0] - xoff;
            point[1] = controls->p[num][1] - yoff;
            gwy_selection_set_object(controls->selection, i, point);
            num++;
        }
    }
    preview(controls);
}

static void
zoom_adjust_peaks(ThresholdControls *controls)
{
    int i, num = 0;
    for (i = 0; i < 4; i++)
    {
        double point[2];
        if (gwy_selection_get_object(controls->selection, i, point))
        {
            gdouble xoff, yoff, multiplier;
            multiplier = 1.0/controls->args->zoom_mode;
            xoff = gwy_data_field_get_xoffset(controls->corr_fft) * multiplier;
            yoff = gwy_data_field_get_yoffset(controls->corr_fft) * multiplier;
            point[0] = controls->p[num][0] - xoff;
            point[1] = controls->p[num][1] - yoff;
            gwy_selection_set_object(controls->selection, i, point);
            num++;
        }
    }
    reFind_Peaks(controls);
}

static void
skew_process(ThresholdControls *controls)
{
    GwyDataField *temp;
    gint i;
    gdouble Trans[6], iTrans[6];
    gdouble p[3], Tp[3], cornX[4], cornY[4], TcornX[4], TcornY[4];
    gdouble lowX, highX, lowY, highY;
    gdouble oxres, oyres, xres, yres;
    gdouble xreal, yreal, xscale, yscale;
    gdouble min, max, hAngle, vAngle;
    g_object_unref(controls->corr_image);
    g_object_unref(controls->corr_fft);
    oxres = gwy_data_field_get_xres(controls->image);
    oyres = gwy_data_field_get_yres(controls->image);
    gwy_data_field_get_min_max(controls->image, &min, &max);
    controls->args->background_fill = min - 0.05 * (max - min);
    hAngle = deg2rad(controls->args->Xskew);
    vAngle = deg2rad(controls->args->Yskew);
    cornX[0] = 0;
    cornX[1] = oxres;
    cornX[2] = oxres;
    cornX[3] = 0;
    cornY[0] = 0;
    cornY[1] = 0;
    cornY[2] = oyres;
    cornY[3] = oyres;
    Trans[0] = 1;
    Trans[1] = tan(vAngle);
    Trans[2] = tan(hAngle);
    Trans[3] = 1;
    Trans[4] = 0;
    Trans[5] = 0;
    for (i = 0; i < 4; i++)
    {
        p[0] = cornX[i];
        p[1] = cornY[i];
        p[2] = 1;
        mult_3matrix(Tp, Trans, p);
        TcornX[i] = Tp[0];
        TcornY[i] = Tp[1];
    }
    lowX = TcornX[0];
    highX = TcornX[0];
    lowY = TcornY[0];
    highY = TcornY[0];
    for (i = 1; i < 4; i++)
    {
        if (TcornX[i] < lowX)
            lowX = TcornX[i];
        if (TcornX[i] > highX)
            highX = TcornX[i];
        if (TcornY[i] < lowY)
            lowY = TcornY[i];
        if (TcornY[i] > highY)
            highY = TcornY[i];
    }
    xres = GWY_ROUND(highX - lowX);
    yres = GWY_ROUND(highY - lowY);
    controls->args->newxres = xres;
    controls->args->newyres = yres;
    xscale = xres/oxres;
    yscale = yres/oyres;
    xreal = gwy_data_field_get_xreal(controls->image) * xscale;
    yreal = gwy_data_field_get_yreal(controls->image) * yscale;
    controls->corr_image = gwy_data_field_new(xres, yres, xreal, yreal, FALSE);
    gwy_data_field_fill(controls->corr_image, controls->args->background_fill);
    temp = gwy_data_field_duplicate(controls->image);
    Trans[0] = 1;
    Trans[1] = tan(vAngle);
    Trans[2] = tan(hAngle);
    Trans[3] = 1;
    Trans[4] = -lowX;
    Trans[5] = -lowY;
    invert_matrix(iTrans, Trans);
    affine(temp, controls->corr_image, iTrans, GWY_INTERPOLATION_BILINEAR,
            controls->args->background_fill);
    g_object_unref(temp);
    gwy_data_field_set_si_unit_xy(controls->corr_image,
            controls->Image_XY_Units);
    gwy_data_field_set_si_unit_z(controls->corr_image,
            controls->Image_Z_Units);
    controls->corr_fft = gwy_data_field_duplicate(controls->corr_image);
    perform_fft(controls->corr_fft, controls->mydata);
}

static void
skew_create_output(GwyContainer *data,
    GwyDataField *dfield, ThresholdControls *controls)
{
    const guchar *title;
    GwyContainer *meta;
    gint id, newid;
    gdouble oxres, oyres, xres, yres;
    gdouble xreal, yreal, xscale, yscale;
    oxres = gwy_data_field_get_xres(controls->image);
    oyres = gwy_data_field_get_yres(controls->image);
    xres = controls->args->newxres;
    yres = controls->args->newyres;
    xscale = xres/oxres;
    yscale = yres/oyres;
    xreal = gwy_data_field_get_xreal(controls->image) * xscale;
    yreal = gwy_data_field_get_yreal(controls->image) * yscale;
    gwy_data_field_set_xreal(dfield, xreal);
    gwy_data_field_set_yreal(dfield, yreal);
    gwy_data_field_set_si_unit_xy(dfield,
        controls->Image_XY_Units);
    gwy_data_field_set_si_unit_z(dfield,
        controls->Image_Z_Units);
    gwy_app_data_browser_get_current(GWY_APP_DATA_FIELD_ID, &id, 0);
    GQuark Qmeta = g_quark_from_string(g_strdup_printf("/%i/meta", id));
    if (gwy_container_contains(data, Qmeta))
        meta = gwy_container_duplicate(gwy_container_get_object(data, Qmeta));
    else
        meta = gwy_container_new();
    title = gwy_container_get_string(data,
            g_quark_try_string(g_strdup_printf("/%i/data/title", id)));
    gwy_container_set_string_by_name(meta, "Source Title", title);
    gwy_container_set_string_by_name(meta, "X Skew (°)",
            (const guchar *)g_strdup_printf("%.5f", controls->args->Xskew));
    gwy_container_set_string_by_name(meta, "Y Skew (°)",
            (const guchar *)g_strdup_printf("%.5f", controls->args->Yskew));
    newid = gwy_app_data_browser_add_data_field(dfield, data, TRUE);
    gwy_container_set_object_by_name(data,
            g_strdup_printf("/%i/meta", newid), meta);
    gwy_app_set_data_field_title(data, newid, _("Skewed"));
    gwy_app_channel_log_add(data, controls->id,
        newid, "proc::skew_lattice", NULL);
}

static void
gwy_tool_level3_render_cell(GtkCellLayout *layout,
            GtkCellRenderer *renderer, GtkTreeModel *model,
            GtkTreeIter *iter, gpointer user_data)
{
    ThresholdControls *controls = (ThresholdControls*)user_data;
    const GwySIValueFormat *vf;
    gchar buf[32];
    gdouble point[2];
    gdouble val;
    guint idx, id;
    id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(layout), "id"));
    gtk_tree_model_get(model, iter, 0, &idx, -1);
    if (id == COLUMN_I)
    {
        g_snprintf(buf, sizeof(buf), "%d", idx + 1);
        g_object_set(renderer, "text", buf, NULL);
        return;
    }
    if (!controls->selection ||
            !gwy_selection_get_object(controls->selection, idx, point))
    {
        g_object_set(renderer, "text", "", NULL);
        return;
    }
    switch (id)
    {
        case COLUMN_X:
            if (controls->args->image_mode == IMAGE_FFT_CORRECTED)
                peak_find(controls, point, idx);
            vf = controls->XY_Format;
            val = controls->p[idx][0];
            break;
        case COLUMN_Y:
            vf = controls->XY_Format;
            val = controls->p[idx][1];
            break;
        case COLUMN_Z:
            vf = controls->Z_Format;
            val = controls->p[idx][2];
            break;
        default:
            g_return_if_reached();
            break;
    }
    if (vf)
        g_snprintf(buf, sizeof(buf), "%.*f",
            vf->precision, val/vf->magnitude);
    else
        g_snprintf(buf, sizeof(buf), "%.3g", val);
    g_object_set(renderer, "text", buf, NULL);
    skew_update_angles(controls);
}

static void
skew_Xadjusted(ThresholdControls *controls)
{
    GtkAdjustment *adj = (GtkAdjustment*)controls->skew_Xadjust;
    controls->args->Xskew = adj->value;
    skew_process(controls);
    reFind_Peaks(controls);
    gchar *s = g_strdup_printf("%0.1f", adj->value);
    gtk_entry_set_text(GTK_ENTRY(controls->hskewtxt), s);
    g_free(s);
}

static void
skew_Yadjusted(ThresholdControls *controls)
{
    GtkAdjustment *adj = (GtkAdjustment*)controls->skew_Yadjust;
    controls->args->Yskew = adj->value;
    skew_process(controls);
    reFind_Peaks(controls);
    gchar *s = g_strdup_printf("%0.1f", adj->value);
    gtk_entry_set_text(GTK_ENTRY(controls->vskewtxt), s);
    g_free(s);
}

static void
skew_update_angles(ThresholdControls *controls)
{
    if (gwy_selection_is_full(controls->selection))
    {
        get_angles(controls);
        gtk_label_set_markup(GTK_LABEL(controls->Angle1),
            g_strdup_printf("<b>Angle 123: </b>%0.1f°",
                                        controls->args->angle1));
        gtk_label_set_markup(GTK_LABEL(controls->Angle2),
            g_strdup_printf("<b>Angle 234: </b>%0.1f°",
                                        controls->args->angle2));
    }
    else
    {
        gtk_label_set_markup(GTK_LABEL(controls->Angle1), "<b>Angle 123:</b>");
        gtk_label_set_markup(GTK_LABEL(controls->Angle2), "<b>Angle 234:</b>");
    }
}

static void
get_angles(ThresholdControls *controls)
{
    gdouble x1, x2, x3, x4, y1, y2, y3, y4;
    gdouble ax, ay, bx, by, a, b, ab;
    x1 = controls->p[0][0];
    y1 = controls->p[0][1];
    x2 = controls->p[1][0];
    y2 = controls->p[1][1];
    x3 = controls->p[2][0];
    y3 = controls->p[2][1];
    x4 = controls->p[3][0];
    y4 = controls->p[3][1];
    ax = x1 - x2;
    ay = y1 - y2;
    bx = x3 - x2;
    by = y3 - y2;
    a = sqrt(ax*ax + ay*ay);
    b = sqrt(bx*bx + by*by);
    ab = ax*bx + ay*by;
    controls->args->angle1 = acos(ab / (a*b)) * 180/PI;
    ax = x4 - x3;
    ay = y4 - y3;
    bx = x2 - x3;
    by = y2 - y3;
    a = sqrt(ax*ax + ay*ay);
    ab = ax*bx + ay*by;
    controls->args->angle2 = acos(ab / (a*b)) * 180/PI;
}

static void
reset_Xskew(ThresholdControls *controls)
{
    gtk_adjustment_set_value((GtkAdjustment*)controls->skew_Xadjust, 0.0);
}

static void
reset_Yskew(ThresholdControls *controls)
{
    gtk_adjustment_set_value((GtkAdjustment*)controls->skew_Yadjust, 0.0);
}

static void
gwy_tool_level3_radius_changed(GwyToolLevel3 *tool)
{
    tool->rpx = gwy_adjustment_get_int(tool->radius);
    guint i;
    GwyNullStore *store = GWY_NULL_STORE(tool->model);
    for (i = 0; i < 4; i++)
        gwy_null_store_row_changed(store, i);
}

static void
perform_fft(GwyDataField *dfield, GwyContainer *data)
{    
    GwyDataField *raout, *ipout;
    raout = gwy_data_field_new_alike(dfield, FALSE);
    ipout = gwy_data_field_new_alike(dfield, FALSE);
    gwy_data_field_2dfft(dfield, NULL, raout, ipout,
                         GWY_WINDOWING_HANN,
                         GWY_TRANSFORM_DIRECTION_FORWARD,
                         GWY_INTERPOLATION_LINEAR, FALSE, 1);
    set_dfield_modulus(raout, ipout, dfield);
    fft_postprocess(dfield);
    gchar *key;
    key = g_strdup_printf("/%i/base/palette", 0);
    gwy_container_set_string_by_name(data, key, g_strdup("Gray"));
    g_free(key);
    key = g_strdup_printf("/%i/base/range-type", 0);
    gwy_container_set_enum_by_name(data, key, GWY_LAYER_BASIC_RANGE_ADAPT);
    g_free(key);
    g_object_unref(raout);
    g_object_unref(ipout);
}

static void
set_dfield_modulus(GwyDataField *re, GwyDataField *im, GwyDataField *target)
{
    const gdouble *datare, *dataim;
    gdouble *data;
    gint xres, yres, i;
    xres = gwy_data_field_get_xres(re);
    yres = gwy_data_field_get_yres(re);
    datare = gwy_data_field_get_data_const(re);
    dataim = gwy_data_field_get_data_const(im);
    data = gwy_data_field_get_data(target);
    for (i = xres*yres; i; i--, datare++, dataim++, data++)
        *data = hypot(*datare, *dataim);
}

static void
fft_postprocess(GwyDataField *dfield)
{
    gint res;
    gdouble r;
    gwy_data_field_2dfft_humanize(dfield);
    GwySIUnit *xyunit;
    xyunit = gwy_data_field_get_si_unit_xy(dfield);
    gwy_si_unit_power(xyunit, -1, xyunit);
    gwy_data_field_set_xreal(dfield, 1.0/gwy_data_field_get_xmeasure(dfield));
    gwy_data_field_set_yreal(dfield, 1.0/gwy_data_field_get_ymeasure(dfield));
    res = gwy_data_field_get_xres(dfield);
    r = res / 2.0;
    gwy_data_field_set_xoffset(dfield, -gwy_data_field_jtor(dfield, r));
    res = gwy_data_field_get_yres(dfield);
    r = res / 2.0;
    gwy_data_field_set_yoffset(dfield, -gwy_data_field_itor(dfield, r));
    gdouble dmin, dmax;
    gwy_data_field_get_min_max(dfield, &dmin, &dmax);
    gwy_data_field_add(dfield, -dmin);
}

static void
selection_changed(ThresholdControls *controls)
{
    guint i;
    GwyNullStore *store = GWY_NULL_STORE(controls->tool->model);
    for (i = 0; i < 4; i++)
        gwy_null_store_row_changed(store, i);
}

static void
skew_do(ThresholdControls *controls)
{
    skew_process(controls);
    skew_create_output(controls->container, controls->corr_image, controls);
    g_object_unref(controls->image);
    g_object_unref(controls->dfield);
    g_object_unref(controls->corr_image);
    g_object_unref(controls->corr_fft);
}

static void
image_mode_changed(GtkToggleButton *button, ThresholdControls *controls)
{
    threshold_load_args(controls);
    controls->args->image_mode =
        gwy_radio_buttons_get_current(controls->image_mode_radios);
    if (controls->args->image_mode != IMAGE_FFT_CORRECTED)
        gwy_vector_layer_set_editable(controls->vlayer, FALSE);
    else
        gwy_vector_layer_set_editable(controls->vlayer, TRUE);
    preview(controls);
}

static void
zoom_mode_changed(GtkToggleButton *button, ThresholdControls *controls)
{
    controls->args->zoom_mode =
        gwy_radio_buttons_get_current(controls->zoom_mode_radios);
    preview(controls);
    zoom_adjust_peaks(controls);
}

static void
threshold_do(ThresholdArgs *args, GwyDataField *dfield)
{
    gdouble lower = MIN(args->lower, args->upper);
    gdouble upper = MAX(args->lower, args->upper);
    gwy_data_field_clamp(dfield, lower, upper);
    gwy_data_field_data_changed(dfield);
}

static void
clear_points(ThresholdControls *controls)
{
    gwy_selection_clear(controls->selection);
    skew_update_angles(controls);
    preview(controls);
}

static const gchar lower0_key[] = "/module/skew_lattice/lower0";
static const gchar lower1_key[] = "/module/skew_lattice/lower1";
static const gchar lower2_key[] = "/module/skew_lattice/lower2";
static const gchar lower3_key[] = "/module/skew_lattice/lower3";
static const gchar upper0_key[] = "/module/skew_lattice/upper0";
static const gchar upper1_key[] = "/module/skew_lattice/upper1";
static const gchar upper2_key[] = "/module/skew_lattice/upper2";
static const gchar upper3_key[] = "/module/skew_lattice/upper3";
static const gchar radius_key[] = "/module/skew_lattice/radius";

static void
threshold_load_args(ThresholdControls *controls)
{
    GwyContainer *settings = gwy_app_settings_get();
    gdouble *lower = &controls->args->lower;
    gdouble *upper = &controls->args->upper;
    switch (controls->args->image_mode)
    {
        case IMAGE_DATA:
            gwy_container_gis_double_by_name(settings, lower0_key, lower);
            gwy_container_gis_double_by_name(settings, upper0_key, upper);
            break;
        case IMAGE_FFT:
            gwy_container_gis_double_by_name(settings, lower1_key, lower);
            gwy_container_gis_double_by_name(settings, upper1_key, upper);
            break;
        case IMAGE_CORRECTED:
            gwy_container_gis_double_by_name(settings, lower2_key, lower);
            gwy_container_gis_double_by_name(settings, upper2_key, upper);
            break;
        case IMAGE_FFT_CORRECTED:
            gwy_container_gis_double_by_name(settings, lower3_key, lower);
            gwy_container_gis_double_by_name(settings, upper3_key, upper);
            break;
    }
    threshold_format_value(controls,
                        GTK_ENTRY(controls->upper), controls->args->upper);
    threshold_format_value(controls,
                        GTK_ENTRY(controls->lower), controls->args->lower);
    gwy_container_gis_int32_by_name(settings,
                        radius_key, (&controls->tool->rpx));
}

static void
threshold_save_args(ThresholdControls *controls)
{
    GwyContainer *settings = gwy_app_settings_get();
    switch (controls->args->image_mode)
    {
        case IMAGE_DATA:
            gwy_container_set_double_by_name(settings,
                                lower0_key, controls->args->lower);
            gwy_container_set_double_by_name(settings,
                                upper0_key, controls->args->upper);
            break;
        case IMAGE_FFT:
            gwy_container_set_double_by_name(settings,
                                lower1_key, controls->args->lower);
            gwy_container_set_double_by_name(settings,
                                upper1_key, controls->args->upper);
            break;
        case IMAGE_CORRECTED:
            gwy_container_set_double_by_name(settings,
                                lower2_key, controls->args->lower);
            gwy_container_set_double_by_name(settings,
                                upper2_key, controls->args->upper);
            break;
        case IMAGE_FFT_CORRECTED:
            gwy_container_set_double_by_name(settings,
                                lower3_key, controls->args->lower);
            gwy_container_set_double_by_name(settings,
                                upper3_key, controls->args->upper);
            break;
    }
    gwy_container_set_int32_by_name(settings, radius_key, controls->tool->rpx);
}

static void
hskew_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->hskewtxt));
    gdouble num = g_strtod(value, NULL);
    if (num != controls->args->Xskew)
    {
        controls->args->Xskew = num;
        gtk_adjustment_set_value((GtkAdjustment*)controls->skew_Xadjust, num);
    }
    else
    {
        gchar *s = g_strdup_printf("%0.1f", controls->args->Xskew);
        gtk_entry_set_text(GTK_ENTRY(controls->hskewtxt), s);
        g_free(s);
    }
}

static void
vskew_changed(ThresholdControls *controls)
{
    const gchar *value = gtk_entry_get_text(GTK_ENTRY(controls->vskewtxt));
    gdouble num = g_strtod(value, NULL);
    if (num != controls->args->Yskew)
    {
        controls->args->Yskew = num;
        gtk_adjustment_set_value((GtkAdjustment*)controls->skew_Yadjust, num);
    }
    else
    {
        gchar *s = g_strdup_printf("%0.1f", controls->args->Yskew);
        gtk_entry_set_text(GTK_ENTRY(controls->vskewtxt), s);
        g_free(s);
    }
}

static gdouble
matrix_det(const gdouble *m)
{
    return m[0]*m[3] - m[1]*m[2];
}

static void
invert_matrix(gdouble *dest, const gdouble *src)
{
    gdouble D = matrix_det(src);
    dest[0] = src[3]/D;
    dest[1] = -src[1]/D;
    dest[2] = -src[2]/D;
    dest[3] = src[0]/D;
    dest[4] = (src[2]*src[5] - src[3]*src[4])/D;
    dest[5] = (src[1]*src[4] - src[0]*src[5])/D;
}

static void
mult_3matrix(gdouble *dest, const gdouble *mat1, const gdouble *mat2)
{
    dest[0] = mat1[0]*mat2[0] + mat1[2]*mat2[1] + mat1[4]*mat2[2];
    dest[1] = mat1[1]*mat2[0] + mat1[3]*mat2[1] + mat1[5]*mat2[2];
    dest[2] = mat2[2];
}

static gdouble
deg2rad(const gdouble deg)
{
    return deg * PI / 180.0;
}

static void
affine(GwyDataField *source, GwyDataField *dest, const gdouble *invtrans,
            GwyInterpolationType interp, gdouble fill_value)
{
    GwyDataField *coeffield;
    gdouble *data, *coeff;
    const gdouble *cdata;
    gint xres, yres, newxres, newyres;
    gint newi, newj, oldi, oldj, i, j, ii, jj, suplen, sf, st;
    gdouble x, y, v;
    gdouble axx, axy, ayx, ayy, bx, by;
    gboolean vset;
    g_return_if_fail(GWY_IS_DATA_FIELD(source));
    g_return_if_fail(GWY_IS_DATA_FIELD(dest));
    g_return_if_fail(invtrans);
    axx = invtrans[0];
    axy = invtrans[1];
    ayx = invtrans[2];
    ayy = invtrans[3];
    bx = invtrans[4];
    by = invtrans[5];
    suplen = gwy_interpolation_get_support_size(interp);
    g_return_if_fail(suplen > 0);
    coeff = g_newa(gdouble, suplen*suplen);
    sf = -((suplen - 1)/2);
    st = suplen/2;
    xres = gwy_data_field_get_xres(source);
    yres = gwy_data_field_get_yres(source);
    newxres = gwy_data_field_get_xres(dest);
    newyres = gwy_data_field_get_yres(dest);
    if (gwy_interpolation_has_interpolating_basis(interp))
        coeffield = g_object_ref(source);
    else
    {
        coeffield = gwy_data_field_duplicate(source);
        gwy_interpolation_resolve_coeffs_2d(xres, yres, xres,
                            gwy_data_field_get_data(coeffield),
                            interp);
    }
    data = gwy_data_field_get_data(dest);
    cdata = gwy_data_field_get_data_const(coeffield);
    bx += 0.5*(axx + axy - 1.0);
    by += 0.5*(ayx + ayy - 1.0);
    for (newi = 0; newi < newxres; newi++)
    {
        for (newj = 0; newj < newyres; newj++)
        {
            x = axx*newi + ayx*newj + bx;
            y = axy*newi + ayy*newj + by;
            vset = FALSE;
            if (y > yres || x > xres || y < 0.0 || x < 0.0) {
                v = fill_value;
                vset = TRUE;
            }
            if (!vset) {
                oldi = (gint)floor(y);
                y -= oldi;
                oldj = (gint)floor(x);
                x -= oldj;
                for (i = sf; i <= st; i++) {
                    ii = (oldi + i + 2*st*yres) % (2*yres);
                    if (G_UNLIKELY(ii >= yres))
                        ii = 2*yres-1 - ii;
                    for (j = sf; j <= st; j++) {
                        jj = (oldj + j + 2*st*xres) % (2*xres);
                        if (G_UNLIKELY(jj >= xres))
                            jj = 2*xres-1 - jj;
                        coeff[(i - sf)*suplen + j - sf] = cdata[ii*xres + jj];
                    }
                }
                v = gwy_interpolation_interpolate_2d(x, y, suplen, coeff,
                                                     interp);
            }
            data[newi + newxres*newj] = v;
        }
    }
    g_object_unref(coeffield);
}
