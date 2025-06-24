#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <locale.h>
#ifdef __APPLE__
#include <xlocale.h>
#endif

#include <yaml.h>

#include "demes.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

/*
 * The structs for the defaults aren't in demes.h because they're not
 * needed once the graph is fully resolved.
 */
struct toplevel_defaults {
    /* defaults.epoch */
    struct demes_epoch epoch;
    /* defaults.deme
     * We don't reuse struct demes_deme from demes.h, because the name and
     * epochs fields aren't allowed in defaults.deme, and in defaults.deme
     * we must consider the possibility of different lengths for the ancestors
     * and proportions lists. Finally, we must store the ancestors as a list
     * of names here, because they're resolved only when the demes are parsed.
     */
    struct toplevel_defaults_deme {
        demes_char_t *description;
        double start_time;
        demes_char_t **ancestors; // list of deme names
        double *proportions;
        size_t n_ancestors;
        size_t n_proportions;
    } deme;
    /* defaults.migration
     * We store the source and dest as deme names here, because the're resolved
     * later when the migrations list is parsed. Furthermore, we need a list
     * of deme names for the 'demes' field, which is not a field of
     * struct demes_migration from demes.h.
     */
    struct toplevel_defaults_migration {
        demes_char_t *source;
        demes_char_t *dest;
        double start_time;
        double end_time;
        double rate;
        // list of deme names for symmetric migration
        demes_char_t **demes;
        size_t n_demes;
    } migration;
    /* defaults.pulse
     * Just as for demes.migration above, the source and dest are stored as
     * deme names here.
     */
    struct toplevel_defaults_pulse {
        demes_char_t **sources;
        demes_char_t *dest;
        double time;
        double *proportions;
        size_t n_sources;
        size_t n_proportions;
    } pulse;
};

struct demelevel_defaults {
    struct demes_epoch epoch;
};

/*
 * Print a possibly-UTF8 message to stderr.
 * TODO:
 *      * Some terminals don't support UTF8. There doesn't appear to be any
 *        sane, standard or portable way to print non-ASCII chars though.
 *      * This is a library, so we shouldn't print anything to the terminal.
 */
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif
static void
errmsg(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

/**
 * @return 1 if a and b are close, 0 otherwise.
 */
static int
is_close(double a, double b) {
    return fabs(b - a) < 1e-9;
}

/*
 * Sets out_codepoint to the current utf8 codepoint in str, and returns the
 * address of the next utf8 codepoint after the current one in str.
 *
 * Stolen from https://github.com/sheredom/utf8.h/blob/master/utf8.h
 */
static void *
utf8codepoint(const void *restrict str, uint32_t *restrict out_codepoint)
{
    const char *s = (const char *)str;

    if (0xf0 == (0xf8 & s[0])) {
        // 4 byte utf8 codepoint
        *out_codepoint = ((0x07 & s[0]) << 18) | ((0x3f & s[1]) << 12) |
                         ((0x3f & s[2]) << 6) | (0x3f & s[3]);
        s += 4;
    } else if (0xe0 == (0xf0 & s[0])) {
        // 3 byte utf8 codepoint
        *out_codepoint =
            ((0x0f & s[0]) << 12) | ((0x3f & s[1]) << 6) | (0x3f & s[2]);
        s += 3;
    } else if (0xc0 == (0xe0 & s[0])) {
        // 2 byte utf8 codepoint
        *out_codepoint = ((0x1f & s[0]) << 6) | (0x3f & s[1]);
        s += 2;
    } else {
        // 1 byte utf8 codepoint otherwise
        *out_codepoint = s[0];
        s += 1;
    }

    return (void *)s;
}

/**
 * Check if the given deme name is valid.
 *
 * Deme names must be python identifiers. Lots of unicode is allowed.
 * https://docs.python.org/3/reference/lexical_analysis.html#identifiers
 *
 * @return 1 if the name is a valid deme name, 0 otherwise.
 */
static int
is_valid_deme_name(demes_char_t *name)
{
    demes_char_t *c = name;
    uint32_t ucp; // unicode code point
    extern int _PyUnicode_IsXidStart(uint32_t ch);
    extern int _PyUnicode_IsXidContinue(uint32_t ch);

    if (c == NULL || *c == '\0') {
        // empty string is not valid
        return 0;
    }

    // leading char
    c = utf8codepoint(c, &ucp);
    if (!_PyUnicode_IsXidStart(ucp) && ucp != '_') {
        return 0;
    }

    while (*c != '\0') {
        // subsequent chars
        c = utf8codepoint(c, &ucp);
        if (!_PyUnicode_IsXidContinue(ucp)) {
            return 0;
        }
    }
    return 1;
}

/* strdup wrapper for demes_char_t (unsigned char) */
static inline demes_char_t *
u8_strdup(const demes_char_t *str)
{
    return (demes_char_t *)strdup((const char *)str);
}

/* strncmp wrapper for demes_char_t (unsigned char) */
static inline int
u8_strcmp(const demes_char_t *s1, const demes_char_t *s2)
{
    return strcmp((const char *)s1, (const char *)s2);
}

/**
 * Copy the list of strings from source into dest.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
u8_strdup_list(demes_char_t ***_dest, demes_char_t **source, size_t n_strings)
{
    demes_char_t **dest;
    int ret = 0;
    int i;

    *_dest = NULL;

    dest = calloc(n_strings, sizeof *dest);
    if (dest == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for (i=0; i<n_strings; i++) {
        if ((dest[i] = u8_strdup(source[i])) == NULL) {
            perror("strdup");
            ret = DEMES_ERR_MEMORY;
            goto err1;
        }
    }

    *_dest = dest;

err1:
    if (ret) {
        for (i=0; i<n_strings; i++) {
            if (dest[i]) {
                free(dest[i]);
            }
        }
    }
err0:
    return ret;
}

/**
 * Append an epoch to the @p deme.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_deme_add_epoch(
    struct demes_deme *deme, double end_time, double start_size, double end_size,
    enum size_function size_function, double selfing_rate, double cloning_rate)
{
    struct demes_epoch *epoch;
    int ret = 0;
    double start_time;

    if (isinf(end_time) || end_time < 0) {
        errmsg("deme %s: epochs[%zd]: must have infinity > end_time >= 0\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (deme->n_epochs == 0) {
        start_time = deme->start_time;
    } else {
        start_time = deme->epochs[deme->n_epochs - 1].end_time;
    }
    if (start_time <= end_time) {
        errmsg("deme %s: epochs[%zd]: must have epoch start_time > end_time\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (isinf(start_size) || start_size <= 0) {
        errmsg("deme %s: epochs[%zd]: must have infinity > start_size > 0\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (isinf(end_size) || end_size <= 0) {
        errmsg("deme %s: epochs[%zd]: must have infinity > end_size > 0\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    assert(size_function == DEMES_SIZE_FUNCTION_CONSTANT ||
            size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL);

    if (deme->n_epochs == 0 && isinf(deme->start_time)
            && size_function != DEMES_SIZE_FUNCTION_CONSTANT) {
        errmsg("deme %s: if deme start_time is infinite, size must be constant "
                "in the first epoch\n", deme->name);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (size_function == DEMES_SIZE_FUNCTION_CONSTANT && start_size != end_size) {
        errmsg("deme %s: size_function is constant, but start_size != end_size\n",
                deme->name);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (selfing_rate > 1.0 || selfing_rate < 0.0) {
        errmsg("deme %s: epochs[%zd]: must have 0 >= selfing_rate >= 0\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (cloning_rate > 1.0 || cloning_rate < 0.0) {
        errmsg("deme %s: epochs[%zd]: must have 0 >= cloning_rate >= 0\n",
                deme->name, deme->n_epochs);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    epoch = &deme->epochs[deme->n_epochs];
    epoch->end_time = end_time;
    epoch->start_size = start_size;
    epoch->end_size = end_size;
    epoch->size_function = size_function;
    epoch->selfing_rate = selfing_rate;
    epoch->cloning_rate = cloning_rate;
    deme->n_epochs++;

err0:
    return ret;
}

/**
 * Get the index of the @p deme in the @p graph's demes list.
 *
 * @return The index, or -1 if not found.
 */
static int
get_deme_id(struct demes_graph *graph, demes_char_t *name)
{
    int i;
    for (i=0; i<graph->n_demes; i++) {
        if (!u8_strcmp(graph->demes[i].name, name)) {
            return i;
        }
    }
    return -1;
}

/**
 * Get the deme in the @p graph with the given @p name.
 *
 * @return The deme, or NULL if not found.
 */
static struct demes_deme *
demes_graph_get_deme(struct demes_graph *graph, demes_char_t *name)
{
    int i = get_deme_id(graph, name);
    if (i >= 0) {
        return graph->demes + i;
    } else {
        return NULL;
    }
}

/*
 * Add a deme to the @p graph.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_add_deme(
    struct demes_graph *graph,
    demes_char_t *name, demes_char_t *description, double start_time,
    demes_char_t **ancestors, size_t n_ancestors,
    double *proportions, size_t n_proportions)
{
    struct demes_deme *deme;
    int ret = 0;

    if (get_deme_id(graph, name) != -1) {
        errmsg("demes[%zd]: %s already present in the graph.\n",
                graph->n_demes, name);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (!is_valid_deme_name(name)) {
        errmsg("demes[%zd]: '%s' is not a valid deme name. "
                "We recommend choosing a name that starts with a letter or "
                "underscore, and is followed by one or more letters, numbers, "
                "or underscores.\n", graph->n_demes, name);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (start_time <= 0) {
        errmsg("deme %s: start_time must be > 0\n", name);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    assert((ancestors && n_ancestors > 0) ||
            (ancestors == NULL && n_ancestors == 0 && n_proportions == 0));
    assert((proportions && n_proportions > 0) ||
            (proportions == NULL && n_proportions == 0));

    if (ancestors) {
        int i;

        if (n_ancestors != n_proportions) {
            if (n_ancestors == 1 && n_proportions == 0) {
                assert(proportions == NULL);
            } else {
                errmsg("deme %s: ancestors list and proportions list don't match\n",
                        name);
                ret = DEMES_ERR_PROPORTIONS;
                goto err0;
            }
        }

        for (i=0; i<n_ancestors; i++) {
            struct demes_deme *ancestor;
            double anc_end_time;
            int j;

            if ((ancestor = demes_graph_get_deme(graph, ancestors[i])) == NULL) {
                errmsg("deme %s: ancestor deme '%s' not found. "
                        "Note: ancestor demes must be specified "
                        "before their children.\n",
                        name, ancestors[0]);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
            anc_end_time = ancestor->epochs[ancestor->n_epochs - 1].end_time;
            if (!(ancestor->start_time > start_time && start_time >= anc_end_time)) {
                errmsg("deme %s: start_time=%lf is outside the interval of "
                        "existance for ancestor %s (%lf, %lf]\n",
                        name, start_time, ancestor->name, ancestor->start_time,
                        anc_end_time);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
            for (j=0; j<i; j++) {
                if (!u8_strcmp(ancestor->name, ancestors[j])) {
                    errmsg("demes %s: duplicate ancestor '%s'\n", name, ancestors[j]);
                    ret = DEMES_ERR_VALUE;
                    goto err0;
                }
            }
        }

        if (n_proportions > 0) {
            double sum = 0;
            for (i=0; i<n_proportions; i++) {
                if (isnan(proportions[i]) || proportions[i] < 0
                        || proportions[i] > 1) {
                    errmsg("deme %s: must have 0 >= p >= 1 for each proportion p\n",
                            name);
                    ret = DEMES_ERR_PROPORTIONS;
                    goto err0;
                }
                sum += proportions[i];
            }
            if (!is_close(sum, 1.0)) {
                errmsg("deme %s: proportions must sum to 1.0\n", name);
                ret = DEMES_ERR_PROPORTIONS;
                goto err0;
            }
        }
    }

    if (!ancestors && !isinf(start_time)) {
        errmsg("deme %s: finite start_time, but no ancestors specified\n", name);
        ret = DEMES_ERR_MISSING_REQUIRED;
        goto err0;
    }

    deme = &graph->demes[graph->n_demes];

    if ((deme->name = u8_strdup(name)) == NULL) {
        perror("strdup");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }
    if (description) {
        if ((deme->description = u8_strdup(description)) == NULL) {
            perror("strdup");
            ret = DEMES_ERR_MEMORY;
            goto err0;
        }
    }

    deme->start_time = start_time;

    if (ancestors) {
        int i;

        assert(n_ancestors > 0);

        deme->n_ancestors = n_ancestors;
        deme->ancestors = calloc(n_ancestors, sizeof *deme->ancestors);
        if (deme->ancestors == NULL) {
            perror("calloc");
            ret = DEMES_ERR_MEMORY;
            goto err0;
        }
        deme->proportions = calloc(n_ancestors, sizeof *deme->proportions);
        if (deme->proportions  == NULL) {
            perror("calloc");
            ret = DEMES_ERR_MEMORY;
            goto err0;
        }

        for (i=0; i<n_ancestors; i++) {
            *(deme->ancestors + i) = demes_graph_get_deme(graph, ancestors[i]);
            if (n_ancestors == 1) {
                *deme->proportions = 1.0;
            } else {
                *(deme->proportions + i) = proportions[i];
            }
        }
    }

    graph->n_demes++;

err0:
    return ret;
}

/*
 * Get the start and end times of the period when demes d1 and d2 overlap.
 */
static void
time_intersection(
    struct demes_deme *d1, struct demes_deme *d2,
    double *time_start, double *time_end)
{
    *time_start = min(d1->start_time, d2->start_time);
    *time_end = max(d1->epochs[d1->n_epochs - 1].end_time,
                    d2->epochs[d2->n_epochs - 1].end_time);
}

/**
 * Add an asymmetric migration to the @p graph.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_add_migration(
    struct demes_graph *graph,
    demes_char_t *source, demes_char_t *dest,
    double start_time, double end_time, double rate)
{
    struct demes_migration *migration, *migrations;
    struct demes_deme *source_deme, *dest_deme;
    double time_hi, time_lo;
    size_t size;
    int ret = 0;

    if ((source_deme = demes_graph_get_deme(graph, source)) == NULL) {
        errmsg("migrations[%zd]: deme '%s' not in graph\n",
                graph->n_migrations, source);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }
    if ((dest_deme = demes_graph_get_deme(graph, dest)) == NULL) {
        errmsg("migrations[%zd]: deme '%s' not in graph\n",
                graph->n_migrations, dest);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }
    if (source_deme == dest_deme) {
        errmsg("migrations[%zd]: source and dest cannot be the same deme\n",
                graph->n_migrations);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    time_intersection(source_deme, dest_deme, &time_hi, &time_lo);
    if (time_hi <= time_lo) {
        errmsg("migrations[%zd]: deme %s and deme %s never coexist\n",
                graph->n_migrations, source, dest);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (isnan(start_time)) {
        start_time = time_hi;
    } else if (start_time > time_hi || start_time <= time_lo) {
        errmsg("migrations[%zd]: must have %lf >= start_time > %lf, "
                "as defined by the time-intersection of deme %s and deme %s\n",
                graph->n_migrations, time_hi, time_lo, source, dest);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (isnan(end_time)) {
        end_time = time_lo;
    } else if (end_time >= time_hi || end_time < time_lo) {
        errmsg("migrations[%zd]: must have %lf > end_time >= %lf, "
                "as defined by the time-intersection of deme %s and deme %s\n",
                graph->n_migrations, time_hi, time_lo, source, dest);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (end_time >= start_time) {
        errmsg("migrations[%zd]: must have start_time > end_time\n",
                graph->n_migrations);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (rate < 0 || rate > 1) {
        errmsg("migrations[%zd]: must have 0 >= rate >= 1\n",
                graph->n_migrations);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    size = (sizeof *migrations) * (graph->n_migrations + 1);
    migrations = realloc(graph->migrations, size);
    if (migrations == NULL) {
        perror("realloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }
    graph->migrations = migrations;
    graph->n_migrations++;
    migration = &migrations[graph->n_migrations - 1];
    migration->source = source_deme;
    migration->dest = dest_deme;
    migration->start_time = start_time;
    migration->end_time = end_time;
    migration->rate = rate;

err0:
    return ret;
}

/**
 * Add a symmetric migration to the @p graph.
 *
 * Symmetric migrations are resolved by creating two asymmetric migrations
 * for each pair of demes participating in the migration.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_add_symmetric_migration(
    struct demes_graph *graph,
    demes_char_t **demes, size_t n_demes,
    double start_time, double end_time, double rate)
{
    int ret = 0;
    int i, j;

    if (n_demes < 2) {
        errmsg("migration[%zd]: must specify a list of two or more deme names\n",
                graph->n_migrations);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    for (i=0; i<n_demes; i++) {
        if (demes_graph_get_deme(graph, demes[i]) == NULL) {
            errmsg("migration[%zd]: deme '%s' not in graph\n",
                    graph->n_migrations, demes[i]);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
    }

    for (i=0; i<n_demes; i++) {
        for (j=0; j<n_demes; j++) {
            if (i != j) {
                if ((ret = demes_graph_add_migration(graph, demes[i], demes[j],
                                start_time, end_time, rate))) {
                    goto err0;
                }
            }
        }
    }

err0:
    return ret;
}

/*
 * Get a slot in the pulses memory for the given pulse time,
 * such that pulses are insertion sorted in time-descending order.
 * We assume that there is contiguous memory allocated for n_pulses,
 * and that &pulses[n_pulses - 1] is a new unused slot.
 */
struct demes_pulse *
get_pulse_pointer(struct demes_pulse *pulses, size_t n_pulses, double time)
{
    struct demes_pulse *p, *p_last;

    p_last = &pulses[n_pulses - 1];
    for (p=pulses; p!=p_last; p++) {
        if (time > p->time) {
            break;
        }
    }
    // Shift pulses after the given time.
    memmove(p + 1, p, (p_last - p)*(sizeof *p));

    return p;
}

/**
 * Append a pulse to the @p graph.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_add_pulse(
    struct demes_graph *graph, demes_char_t **sources, size_t n_sources,
    demes_char_t *dest, double time, double *proportions, size_t n_proportions)
{
    struct demes_pulse *pulses, *pulse;
    struct demes_deme *source_deme, *dest_deme;
    double time_hi, time_lo;
    double sum;
    size_t size;
    int ret = 0;
    int i, j;

    if (n_sources != n_proportions) {
        errmsg("pulses[%zd]: sources list and proportions list don't match\n",
                graph->n_pulses);
        ret = DEMES_ERR_PROPORTIONS;
        goto err0;
    }

    assert(sources && n_sources > 0);
    assert(proportions && n_proportions > 0);

    if ((dest_deme = demes_graph_get_deme(graph, dest)) == NULL) {
        errmsg("pulses[%zd]: deme '%s' not in graph\n",
                graph->n_pulses, dest);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    for (i=0; i<n_sources; i++) {
        if ((source_deme = demes_graph_get_deme(graph, sources[i])) == NULL) {
            errmsg("pulses[%zd]: source deme '%s' not in graph\n",
                    graph->n_pulses, sources[i]);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
        if (source_deme == dest_deme) {
            errmsg("pulses[%zd]: source and dest cannot be the same deme\n",
                    graph->n_pulses);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        time_intersection(source_deme, dest_deme, &time_hi, &time_lo);
        if (time_hi <= time_lo) {
            errmsg("pulses[%zd]: deme %s and deme %s never coexist\n",
                    graph->n_pulses, sources[i], dest);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if (time > time_hi || time < time_lo) {
            errmsg("pulses[%zd]: must have %lf >= time >= %lf, "
                    "as defined by the time-intersection of deme %s and deme %s\n",
                    graph->n_pulses, time_hi, time_lo, sources[i], dest);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
        if (time == dest_deme->epochs[dest_deme->n_epochs - 1].end_time) {
            errmsg("pulses[%zd]: invalid pulse at time %lf, which is dest "
                    "deme %s's end_time\n",
                    graph->n_pulses, time, dest);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
        if (time == source_deme->start_time) {
            errmsg("pulses[%zd]: invalid pulse at time %lf, which is source "
                    "deme %s's start_time\n",
                    graph->n_pulses, time, sources[i]);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        for (j=0; j<i; j++) {
            if (!u8_strcmp(sources[i], sources[j])) {
                errmsg("pulses[%zd]: duplicate source '%s'\n",
                        graph->n_pulses, sources[i]);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
        }
    }

    /*
     * Check sum of proportions <= 1.
     */
    sum = 0;
    for (i=0; i<n_proportions; i++) {
        if (isnan(proportions[i]) || proportions[i] < 0 || proportions[i] > 1) {
            errmsg("pulses[%zd]: must have 0 >= p >= 1 for each proportion p\n",
                    graph->n_pulses);
            ret = DEMES_ERR_PROPORTIONS;
            goto err0;
        }
        sum += proportions[i];
    }
    if (sum > 1.0) {
        errmsg("pulses[%zd]: sum of pulse proportions > 1 for dest deme %s "
                "at time %lf\n",
                graph->n_pulses, dest, time);
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    size = (sizeof *pulses) * (graph->n_pulses + 1);
    pulses = realloc(graph->pulses, size);
    if (pulses == NULL) {
        perror("realloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }
    graph->pulses = pulses;
    graph->n_pulses++;

    // Insert in time-descending order.
    pulse = get_pulse_pointer(pulses, graph->n_pulses, time);
    pulse->sources = NULL;
    pulse->dest = dest_deme;
    pulse->time = time;
    pulse->proportions = NULL;
    pulse->n_sources = n_sources;

    pulse->sources = calloc(n_sources, sizeof *pulse->sources);
    if (pulse->sources == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }
    pulse->proportions = calloc(n_proportions, sizeof *pulse->proportions);
    if (pulse->proportions  == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for (i=0; i<n_sources; i++) {
        *(pulse->sources + i) = demes_graph_get_deme(graph, sources[i]);
        *(pulse->proportions + i) = proportions[i];
    }

err0:
    return ret;
}

/**
 * Free the memory used by the @p graph.
 */
void
demes_graph_free(struct demes_graph *graph)
{
    int i;

    if (graph == NULL)
        return;

    if (graph->description) {
        free(graph->description);
    }
    if (graph->time_units) {
        free(graph->time_units);
    }
    if (graph->doi) {
        for (i=0; i<graph->n_dois; i++) {
            if (graph->doi[i]) {
                free(graph->doi[i]);
            }
        }
        free(graph->doi);
    }
    if (graph->demes) {
        for (i=0; i<graph->n_demes; i++) {
            struct demes_deme *deme = &graph->demes[i];
            if (deme->name) {
                free(deme->name);
            }
            if (deme->description) {
                free(deme->description);
            }
            if (deme->ancestors) {
                free(deme->ancestors);
            }
            if (deme->proportions) {
                free(deme->proportions);
            }
            if (deme->epochs) {
                free(deme->epochs);
            }
        }
        free(graph->demes);
    }
    if (graph->migrations) {
        free(graph->migrations);
    }
    if (graph->pulses) {
        for (i=0; i<graph->n_pulses; i++) {
            struct demes_pulse *pulse = &graph->pulses[i];
            if (pulse->sources) {
                free(pulse->sources);
            }
            if (pulse->proportions) {
                free(pulse->proportions);
            }
        }
        free(graph->pulses);
    }

    free(graph);
}

/**
 * Allocate and initialise a new graph.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_init(
    struct demes_graph **_graph,
    demes_char_t *description,
    demes_char_t *time_units,
    double generation_time,
    demes_char_t **doi,
    size_t n_dois)
{
    struct demes_graph *graph;
    int ret = 0;

    *_graph = NULL;

    graph = calloc(1, sizeof *graph);
    if (graph == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    if (time_units) {
        graph->time_units = u8_strdup(time_units);
        if (graph->time_units == NULL) {
            perror("strdup");
            ret = DEMES_ERR_MEMORY;
            goto err0;
        }
    } else {
        errmsg("required field 'time_units' not found\n");
        ret = DEMES_ERR_MISSING_REQUIRED;
        goto err0;
    }

    graph->generation_time = 1;
    if (!isnan(generation_time)) {
        if (isinf(generation_time) || generation_time <= 0) {
            errmsg("must have 0 > generation_time > infinity\n");
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
        graph->generation_time = generation_time;
    } else if (u8_strcmp(time_units, (demes_char_t *)"generations")) {
        errmsg("generation_time is a required field when time_units are not generations\n");
        ret = DEMES_ERR_MISSING_REQUIRED;
        goto err0;
    }
    if (!u8_strcmp(time_units, (demes_char_t *)"generations")
            && graph->generation_time != 1) {
        errmsg("time_units are generations but generation_time != 1\n");
        ret = DEMES_ERR_VALUE;
        goto err0;
    }

    if (description) {
        graph->description = u8_strdup(description);
        if (graph->description == NULL) {
            perror("strdup");
            ret = DEMES_ERR_MEMORY;
            goto err0;
        }
    }

    if (doi && n_dois > 0) {
        if ((ret = u8_strdup_list(&graph->doi, doi, n_dois))) {
            goto err0;
        }
    }
    graph->n_dois = n_dois;

    *_graph = graph;

err0:
    if (ret) {
        demes_graph_free(graph);
    }
    return ret;
}

static int
value_is_null(yaml_node_t *value)
{
    if (value->type == YAML_SCALAR_NODE
            && !u8_strcmp(value->data.scalar.value, (demes_char_t*)"null")) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * Get the value node of a (key, value) pair under the given mapping node,
 * which matches the given key_name.
 *
 * @return The yaml_node_t* pointing to the value node.
 */
static yaml_node_t *
get_value(yaml_document_t *document, yaml_node_t *node, char *key_name)
{
    yaml_node_t *value = NULL;
    size_t n_pairs;
    int i;

    assert(node->type == YAML_MAPPING_NODE);

    n_pairs = node->data.mapping.pairs.top - node->data.mapping.pairs.start;

    for (i=0; i<n_pairs; i++) {
        yaml_node_pair_t *pair = node->data.mapping.pairs.start + i;
        yaml_node_t *key = yaml_document_get_node(document, pair->key);

        assert(key != NULL);
        if (key->type == YAML_SCALAR_NODE) {
            if (!u8_strcmp((demes_char_t *)key_name, key->data.scalar.value)) {
                value = yaml_document_get_node(document, pair->value);
                break;
            }
        }
    }

    return value;
}

/**
 * Get the value string of a (key, value) pair under the given mapping node,
 * which matches the given key_name.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
get_string(
    yaml_document_t *document, yaml_node_t *node, char *key_name,
    demes_char_t **string)
{
    yaml_node_t *value;
    int ret = 0;

    *string = NULL;
    value = get_value(document, node, key_name);
    if (value) {
        if (value->type != YAML_SCALAR_NODE || value_is_null(value)) {
            errmsg("line %ld: %s: expected a string\n",
                    value->start_mark.line, key_name);
            ret = DEMES_ERR_TYPE;
        } else {
            *string = value->data.scalar.value;
        }
    }
    return ret;
}

/**
 * Get a list of string values of a (key, value) pair under the given mapping node,
 * which matches the given key_name.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
get_string_list(
    yaml_document_t *document, yaml_node_t *node, char *key_name,
    demes_char_t ***_strings, size_t *_length)
{
    yaml_node_t *list;
    demes_char_t **strings = NULL;
    size_t length = 0;
    int ret = 0;
    int i;

    *_strings = NULL;
    *_length = 0;

    list = get_value(document, node, key_name);
    if (list) {
        if (list->type != YAML_SEQUENCE_NODE) {
            errmsg("line %ld: %s: expected a list\n",
                    list->start_mark.line, key_name);
            ret = DEMES_ERR_TYPE;
            goto err0;
        }

        length = list->data.sequence.items.top - list->data.sequence.items.start;
        if (length > 0) {
            strings = calloc(length, sizeof *strings);
            if (strings == NULL) {
                perror("calloc");
                ret = DEMES_ERR_MEMORY;
                goto err0;
            }
        }

        for (i=0; i<length; i++) {
            yaml_node_t *elm = yaml_document_get_node(
                    document, list->data.sequence.items.start[i]);
            assert(elm != NULL);
            if (elm->type != YAML_SCALAR_NODE || value_is_null(elm)) {
                errmsg("line %ld: %s: expected a list of strings\n",
                        elm->start_mark.line, key_name);
                ret = DEMES_ERR_TYPE;
                goto err1;
            }
            strings[i] = elm->data.scalar.value;
        }
    }

    *_strings = strings;
    *_length = length;

err1:
    if (ret && strings) {
        free(strings);
    }
err0:
    return ret;
}

/**
 * Convert a string to a number.
 *
 * @return The number, or NAN if the string couldn't be converted to a number.
 */
static double
parse_number(char *string)
{
    char *endptr;
    double number;
    int saved_errno = errno;
    
    number = strtod(string, &endptr);
    /* Restore errno. We don't care about underflow or overflow,
     * and its nice to leave no side effects.
     */
    errno = saved_errno;

    if (*endptr != '\0') {
        int sign = 1;
        if (string && string[0] == '-') {
            sign = -1;
            string++;
        } else if (string && string[0] == '+') {
            string++;
        }
        if (!strcmp(string, ".inf") || !strcmp(string, ".Inf") ||
                !strcmp(string, ".INF")) {
            number = sign * INFINITY;
        } else {
            /* No number to convert or trailing garbage.
             * The special value "null" will also be interpreted as NAN.
             */
            number = NAN;
        }
    }

    return number;
}

/**
 * Get the numeric value of a (key, value) pair under the given mapping node,
 * which matches the given key_name.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
get_number(yaml_document_t *document, yaml_node_t *node, char *key_name, double *number)
{
    yaml_node_t *value;
    int ret = 0;

    *number = NAN;
    value = get_value(document, node, key_name);
    if (value) {
        if (value->type == YAML_SCALAR_NODE) {
            *number = parse_number((char *)value->data.scalar.value);
        }
        if (isnan(*number)){
            errmsg("line %ld: %s: expected a number\n",
                    value->start_mark.line, key_name);
            ret = DEMES_ERR_TYPE;
        }
    }

    return ret;
}

/**
 * Get a list of number values of a (key, value) pair under the given mapping node,
 * which matches the given key_name.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
get_number_list(
    yaml_document_t *document, yaml_node_t *node, char *key_name,
    double **_numbers, size_t *_length)
{
    yaml_node_t *list;
    double *numbers = NULL;
    size_t length = 0;
    int ret = 0;
    int i;

    *_numbers = NULL;
    *_length = 0;

    list = get_value(document, node, key_name);
    if (list) {
        if (list->type != YAML_SEQUENCE_NODE) {
            errmsg("line %ld: %s: expected a list\n",
                    list->start_mark.line, key_name);
            ret = DEMES_ERR_TYPE;
            goto err0;
        }

        length = list->data.sequence.items.top - list->data.sequence.items.start;
        if (length > 0) {
            numbers = calloc(length, sizeof *numbers);
            if (numbers == NULL) {
                perror("calloc");
                ret = DEMES_ERR_MEMORY;
                goto err0;
            }
        }

        for (i=0; i<length; i++) {
            yaml_node_t *elm = yaml_document_get_node(
                    document, list->data.sequence.items.start[i]);
            assert(elm != NULL);
            if (elm->type != YAML_SCALAR_NODE || value_is_null(elm)) {
                errmsg("line %ld: %s: expected a list of numbers\n",
                        elm->start_mark.line, key_name);
                ret = DEMES_ERR_TYPE;
                goto err1;
            }
            numbers[i] = parse_number((char *)elm->data.scalar.value);
        }
    }

    *_numbers = numbers;
    *_length = length;

err1:
    if (ret && numbers) {
        free(numbers);
    }
err0:
    return ret;
}

/**
 * Check if the keys under a mapping node are all in the allowed list.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static inline int
check_allowed_mapping(
    yaml_document_t *document, yaml_node_t *node, demes_char_t **allowed, size_t len)
{
    size_t n_pairs;
    int i, j;
    int ret;

    if (node->type != YAML_MAPPING_NODE) {
        errmsg("line %ld: expected key/value pairs.\n", node->start_mark.line);
        ret = DEMES_ERR_TYPE;
        goto err0;
    }

    n_pairs = node->data.mapping.pairs.top - node->data.mapping.pairs.start;

    for (i=0; i<n_pairs; i++) {
        yaml_node_pair_t *pair = node->data.mapping.pairs.start + i;
        yaml_node_t *key = yaml_document_get_node(document, pair->key);
        int key_allowed = 0;

        assert(key != NULL);
        if (key->type != YAML_SCALAR_NODE) {
            errmsg("line %ld: expected key to be a scalar\n", key->start_mark.line);
            ret = DEMES_ERR_KEY;
            goto err0;
        }

        for (j=0; j<len; j++) {
            if (!u8_strcmp(allowed[j], key->data.scalar.value)) {
                key_allowed = 1;
                break;
            }
        }
        if (!key_allowed) {
            errmsg("line %ld: unknown key: \"%s\"\n",
                    key->start_mark.line, key->data.scalar.value);
            ret = DEMES_ERR_KEY;
            goto err0;
        }
    }

    ret = 0;
err0:
    return ret;
}

static int
check_allowed_toplevel(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "description", "time_units", "generation_time", "doi", "demes",
        "migrations", "pulses", "defaults", "metadata"
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_toplevel_defaults(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {"deme", "epoch", "migration", "pulse"};
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_toplevel_defaults_deme(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "description", "start_time", "ancestors", "proportions",
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_demelevel_defaults(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {"epoch"};
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_demeslevel(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "name", "description", "start_time", "ancestors", "proportions",
        "epochs", "defaults"
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_epochlevel(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "end_time", "start_size", "end_size", "size_function",
        "selfing_rate", "cloning_rate",
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_migration(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "demes", "source", "dest", "start_time", "end_time", "rate",
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

static int
check_allowed_pulse(yaml_document_t *document, yaml_node_t *node)
{
    static char *allowed[] = {
        "sources", "dest", "time", "proportions",
    };
    size_t len = sizeof(allowed) / sizeof(allowed[0]);
    return check_allowed_mapping(document, node, (demes_char_t **)allowed, len);
}

/**
 * Parse the demes.epochs node (@p epochs) in the yaml @p document,
 * populating the @p deme->epochs as we go.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_epochs(
    yaml_document_t *document,
    struct demes_graph *graph, struct demes_deme *deme,
    struct demes_epoch *epoch_defaults,
    yaml_node_t *epochs)
{
    size_t n_epochs = 1;
    int ret = 0;
    int i;

    if (epochs) {
        if (epochs->type != YAML_SEQUENCE_NODE) {
            errmsg("line %ld: deme %s: epochs must be a list of epoch objects\n",
                    epochs->start_mark.line, deme->name);
            ret = DEMES_ERR_TYPE;
            goto err0;
        }
        n_epochs = epochs->data.sequence.items.top - epochs->data.sequence.items.start;
        if (n_epochs == 0) {
            /* We could actually allow this case if there were a valid
             * defaults.epoch, but then explicitly specifying an empty list
             * is unnecessary. This is quite weird, so just raise an error.
             */
            errmsg("line %ld: deme %s: epochs must be a list of epoch objects\n",
                    epochs->start_mark.line, deme->name);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
    }

    deme->epochs = calloc(n_epochs, sizeof *deme->epochs);
    if (deme->epochs == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for (i=0; i<n_epochs; i++) {
        double end_time, start_size, end_size, selfing_rate, cloning_rate;
        enum size_function size_function;

        end_time = epoch_defaults->end_time;
        start_size = epoch_defaults->start_size;
        end_size = epoch_defaults->end_size;
        size_function = epoch_defaults->size_function;
        selfing_rate = epoch_defaults->selfing_rate;
        cloning_rate = epoch_defaults->cloning_rate;

        if (epochs) {
            demes_char_t *size_function_str;
            yaml_node_t *epoch = yaml_document_get_node(
                    document, epochs->data.sequence.items.start[i]);

            assert(epoch != NULL);

            if (epoch->type != YAML_MAPPING_NODE) {
                errmsg("line %ld: deme %s: epochs[%d]: expected an epoch object\n",
                        epochs->start_mark.line, deme->name, i);
                ret = DEMES_ERR_TYPE;
                goto err0;
            }
            if ((ret = check_allowed_epochlevel(document, epoch))) {
                goto err0;
            }

            if (get_value(document, epoch, "end_time")) {
                if ((ret = get_number(document, epoch, "end_time", &end_time))) {
                    goto err0;
                }
            }
            if (get_value(document, epoch, "start_size")) {
                if ((ret = get_number(document, epoch, "start_size", &start_size))) {
                    goto err0;
                }
            }
            if (get_value(document, epoch, "end_size")) {
                if ((ret = get_number(document, epoch, "end_size", &end_size))) {
                    goto err0;
                }
            }
            if ((ret = get_string(document, epoch, "size_function", &size_function_str))) {
                goto err0;
            }
            if (size_function_str) {
                if (!u8_strcmp(size_function_str, (demes_char_t *)"constant")) {
                    size_function = DEMES_SIZE_FUNCTION_CONSTANT;
                } else if (!u8_strcmp(size_function_str, (demes_char_t *)"exponential")) {
                    size_function = DEMES_SIZE_FUNCTION_EXPONENTIAL;
                } else {
                    errmsg("deme %s: epochs[%d]: size_function must be constant "
                            "or exponential\n", deme->name, i);
                    ret = DEMES_ERR_VALUE;
                    goto err0;
                }
            }
            if (get_value(document, epoch, "selfing_rate")) {
                if ((ret = get_number(document, epoch, "selfing_rate", &selfing_rate))) {
                    goto err0;
                }
            }
            if (get_value(document, epoch, "cloning_rate")) {
                if ((ret = get_number(document, epoch, "cloning_rate", &cloning_rate))) {
                    goto err0;
                }
            }
        }

        if (isnan(end_time)) {
            if (i == n_epochs - 1) {
                /* end time may be omitted in the final epoch */
                end_time = 0;
            } else {
                errmsg("deme %s: epochs[%d]: required field 'end_time' not found\n",
                        deme->name, i);
                ret = DEMES_ERR_MISSING_REQUIRED;
                goto err0;
            }
        }

        if (i == 0) {
            /* first epoch has special rules */
            if (isnan(start_size) && isnan(end_size)) {
                if (epochs) {
                    errmsg("deme %s: first epoch must have a start_size or end_size\n",
                            deme->name);
                } else {
                    errmsg("deme %s: no epochs defined\n", deme->name);
                }
                ret = DEMES_ERR_MISSING_REQUIRED;
                goto err0;
            }
            if (isnan(start_size)) {
                start_size = end_size;
            }
            if (isnan(end_size)) {
                end_size = start_size;
            }
        } else {
            if (isnan(start_size)) {
                start_size = deme->epochs[i - 1].end_size;
            }
            if (isnan(end_size)) {
                end_size = start_size;
            }
        }

        if (size_function == -1) {
            if (start_size == end_size) {
                size_function = DEMES_SIZE_FUNCTION_CONSTANT;
            } else {
                size_function = DEMES_SIZE_FUNCTION_EXPONENTIAL;
            }
        }

        if (isnan(selfing_rate)) {
            selfing_rate = 0.0;
        }
        if (isnan(cloning_rate)) {
            cloning_rate = 0.0;
        }

        if ((ret = demes_deme_add_epoch(
                deme, end_time, start_size, end_size, size_function,
                selfing_rate, cloning_rate))) {
            goto err0;
        }

    }

err0:
    return ret;
}

/**
 * Parse the defaults.epoch node in the yaml @p document,
 * populating the @p epoch.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_defaults_epoch(
    yaml_document_t *document,
    yaml_node_t *defaults, struct demes_epoch *epoch)
{
    yaml_node_t *epoch_node;
    int ret = 0;

    if ((epoch_node = get_value(document, defaults, "epoch"))) {
        demes_char_t *size_function_str;

        if ((ret = check_allowed_epochlevel(document, epoch_node))) {
            goto err0;
        }

        if ((ret = get_number(document, epoch_node, "end_time",
                        &epoch->end_time))) {
            goto err0;
        }
        if (!isnan(epoch->end_time) &&
                (isinf(epoch->end_time) || epoch->end_time < 0)) {
            yaml_node_t *key = get_value(document, epoch_node, "end_time");
            errmsg("line %ld: must have infinity > end_time >= 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_number(document, epoch_node, "start_size",
                        &epoch->start_size))) {
            goto err0;
        }
        if (!isnan(epoch->start_size) &&
                (isinf(epoch->start_size) || epoch->start_size <= 0)) {
            yaml_node_t *key = get_value(document, epoch_node, "start_size");
            errmsg("line %ld: must have infinity > start_size > 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_number(document, epoch_node, "end_size",
                        &epoch->end_size))) {
            goto err0;
        }
        if (!isnan(epoch->end_size) &&
                (isinf(epoch->end_size) || epoch->end_size <= 0)) {
            yaml_node_t *key = get_value(document, epoch_node, "end_size");
            errmsg("line %ld: must have infinity > end_size > 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_string(document, epoch_node, "size_function",
                        &size_function_str))) {
            goto err0;
        }
        if (size_function_str) {
            if (!u8_strcmp(size_function_str, (demes_char_t *)"constant")) {
                epoch->size_function = DEMES_SIZE_FUNCTION_CONSTANT;
            }
            else if (!u8_strcmp(size_function_str, (demes_char_t *)"exponential")) {
                epoch->size_function = DEMES_SIZE_FUNCTION_EXPONENTIAL;
            } else {
                yaml_node_t *key = get_value(document, epoch_node, "size_function");
                errmsg("line %ld: size_function must be constant or exponential\n",
                        key->start_mark.line);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
        }

        if ((ret = get_number(document, epoch_node, "selfing_rate",
                        &epoch->selfing_rate))) {
            goto err0;
        }
        if (!isnan(epoch->selfing_rate) &&
                (epoch->selfing_rate > 1.0 || epoch->selfing_rate < 0.0)) {
            yaml_node_t *key = get_value(document, epoch_node, "selfing_rate");
            errmsg("line %ld: must have 0 >= selfing_rate >= 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_number(document, epoch_node, "cloning_rate",
                        &epoch->cloning_rate))) {
            goto err0;
        }
        if (!isnan(epoch->cloning_rate) &&
                (epoch->cloning_rate > 1.0 || epoch->cloning_rate < 0.0)) {
            yaml_node_t *key = get_value(document, epoch_node, "cloning_rate");
            errmsg("line %ld: must have 0 >= cloning_rate >= 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
    }

err0:
    return ret;
}

/**
 * Parse the defaults.deme @p node in the yaml @p document,
 * populating the @p deme.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_defaults_deme(
    yaml_document_t *document,
    yaml_node_t *defaults,
    struct toplevel_defaults_deme *deme)
{
    yaml_node_t *deme_node;
    int ret = 0;

    deme_node = get_value(document, defaults, "deme");
    if (deme_node) {

        if ((ret = check_allowed_toplevel_defaults_deme(document, deme_node))) {
            goto err0;
        }

        if ((ret = get_string(document, deme_node, "description",
                        &deme->description))) {
            goto err0;
        }

        if ((ret = get_number(document, deme_node, "start_time",
                        &deme->start_time))) {
            goto err0;
        }
        if (!isnan(deme->start_time) && deme->start_time <= 0) {
            yaml_node_t *key = get_value(document, deme_node, "start_time");
            errmsg("line %ld: must have start_time > 0\n", key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_string_list(document, deme_node, "ancestors",
                        &deme->ancestors, &deme->n_ancestors))) {
            goto err0;
        }

        if ((ret = get_number_list(document, deme_node, "proportions",
                        &deme->proportions, &deme->n_proportions))) {
            goto err0;
        }
        if (deme->proportions) {
            int i;
            double sum = 0;
            for (i=0; i<deme->n_proportions; i++) {
                if (isnan(deme->proportions[i]) || deme->proportions[i] < 0
                        || deme->proportions[i] > 1) {
                    yaml_node_t *key = get_value(document, deme_node, "proportions");
                    errmsg("line %ld: must have 0 >= p >= 1 for each proportion p\n",
                            key->start_mark.line);
                    ret = DEMES_ERR_PROPORTIONS;
                    goto err0;
                }
                sum += deme->proportions[i];
            }
            if (!is_close(sum, 1.0)) {
                yaml_node_t *key = get_value(document, deme_node, "proportions");
                errmsg("line %ld: proportions must sum to 1.0\n", key->start_mark.line);
                ret = DEMES_ERR_PROPORTIONS;
                goto err0;
            }
        }
    }

err0:
    return ret;
}

/**
 * Parse the defaults.migration node in the yaml @p document.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_defaults_migration(
    yaml_document_t *document,
    yaml_node_t *defaults,
    struct toplevel_defaults_migration *migration)
{
    yaml_node_t *migration_node;
    int ret = 0;


    migration_node = get_value(document, defaults, "migration");
    if (migration_node) {
        if ((ret = check_allowed_migration(document, migration_node))) {
            goto err0;
        }

        if ((ret = get_string_list(document, migration_node, "demes",
                        &migration->demes, &migration->n_demes))) {
            goto err0;
        }
        if ((ret = get_string(document, migration_node, "source",
                        &migration->source))) {
            goto err0;
        }
        if ((ret = get_string(document, migration_node, "dest",
                        &migration->dest))) {
            goto err0;
        }

        if ((ret = get_number(document, migration_node, "start_time",
                        &migration->start_time))) {
            goto err0;
        }
        if (!isnan(migration->start_time) && migration->start_time <= 0) {
            yaml_node_t *key = get_value(document, migration_node, "start_time");
            errmsg("line %ld: must have start_time > 0\n", key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_number(document, migration_node, "end_time",
                        &migration->end_time))) {
            goto err0;
        }
        if (!isnan(migration->end_time) && (isinf(migration->end_time)
                || migration->end_time < 0)) {
            yaml_node_t *key = get_value(document, migration_node, "end_time");
            errmsg("line %ld: must have infininty > end_time >= 0\n",
                    key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if ((ret = get_number(document, migration_node, "rate",
                        &migration->rate))) {
            goto err0;
        }
        if (!isnan(migration->rate) && (migration->rate < 0
                || migration->rate > 1)) {
            yaml_node_t *key = get_value(document, migration_node, "rate");
            errmsg("line %ld: must have 0 <= rate <= 1\n", key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }
    }

err0:
    return ret;
}

/**
 * Parse the defaults.pulse node in the yaml @p document,
 * populating @p pulse as we go.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_defaults_pulse(
    yaml_document_t *document,
    yaml_node_t *defaults,
    struct toplevel_defaults_pulse *pulse)
{
    yaml_node_t *pulse_node;
    int ret = 0;

    pulse_node = get_value(document, defaults, "pulse");
    if (pulse_node) {
        if ((ret = check_allowed_pulse(document, pulse_node))) {
            goto err0;
        }

        if (get_value(document, pulse_node, "sources")) {
            if ((ret = get_string_list(document, pulse_node, "sources",
                            &pulse->sources, &pulse->n_sources))) {
                goto err0;
            }
            if (pulse->sources == NULL) {
                yaml_node_t *key = get_value(document, pulse_node, "sources");
                errmsg("line %ld: sources must be a non-empty list of deme names\n",
                        key->start_mark.line);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
        }
        if ((ret = get_string(document, pulse_node, "dest", &pulse->dest))) {
            goto err0;
        }

        if ((ret = get_number(document, pulse_node, "time", &pulse->time))) {
            goto err0;
        }
        if (!isnan(pulse->time) && (isinf(pulse->time) || pulse->time <= 0)) {
            yaml_node_t *key = get_value(document, pulse_node, "time");
            errmsg("line %ld: must have infinity > time >= 0\n", key->start_mark.line);
            ret = DEMES_ERR_VALUE;
            goto err0;
        }

        if (get_value(document, pulse_node, "proportions")) {
            if ((ret = get_number_list(document, pulse_node, "proportions",
                            &pulse->proportions, &pulse->n_proportions))) {
                goto err0;
            }
            if (pulse->proportions == NULL) {
                yaml_node_t *key = get_value(document, pulse_node, "proportions");
                errmsg("line %ld: proportions must be a non-empty list of numbers\n",
                        key->start_mark.line);
                ret = DEMES_ERR_VALUE;
                goto err0;
            }
        }
        if (pulse->proportions != NULL) {
            int i;
            double sum = 0;
            for (i=0; i<pulse->n_proportions; i++) {
                if (isnan(pulse->proportions[i]) || pulse->proportions[i] < 0
                        || pulse->proportions[i] > 1) {
                    yaml_node_t *key = get_value(document, pulse_node, "proportions");
                    errmsg("line %ld: must have 0 <= p <= 1 for all proportions p\n",
                            key->start_mark.line);
                    ret = DEMES_ERR_PROPORTIONS;
                    goto err0;
                }
                sum += pulse->proportions[i];
            }
            if (sum > 1) {
                    yaml_node_t *key = get_value(document, pulse_node, "proportions");
                    errmsg("line %ld: proportions must sum to <= 1\n",
                            key->start_mark.line);
                    ret = DEMES_ERR_PROPORTIONS;
                    goto err0;
            }
        }
    }

err0:
    return ret;
}

static void
init_defaults_epoch(struct demes_epoch *epoch)
{
    epoch->end_time = NAN;
    epoch->start_size = NAN;
    epoch->end_size = NAN;
    epoch->size_function = -1;
    epoch->selfing_rate = NAN;
    epoch->cloning_rate = NAN;
}

static void
init_demelevel_defaults(struct demelevel_defaults *defaults)
{
    init_defaults_epoch(&defaults->epoch);
}

static void
init_toplevel_defaults(struct toplevel_defaults *defaults)
{
    /* We use invalid values to signal that no default value was provided. */
    init_defaults_epoch(&defaults->epoch);

    defaults->deme.description = NULL;
    defaults->deme.start_time = NAN;
    defaults->deme.ancestors = NULL;
    defaults->deme.proportions = NULL;
    defaults->deme.n_ancestors = 0;
    defaults->deme.n_proportions = 0;

    defaults->migration.demes = NULL;
    defaults->migration.n_demes = 0;
    defaults->migration.source = NULL;
    defaults->migration.dest = NULL;
    defaults->migration.start_time = NAN;
    defaults->migration.end_time = NAN;
    defaults->migration.rate = NAN;

    defaults->pulse.sources = NULL;
    defaults->pulse.dest = NULL;
    defaults->pulse.time = NAN;
    defaults->pulse.proportions = NULL;
    defaults->pulse.n_sources = 0;
    defaults->pulse.n_proportions = 0;
}

/**
 * Parse a demes.defaults node in the yaml @p document,
 * merging in the toplevel defaults.deme values as appropriate.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_demelevel_defaults(
    yaml_document_t *document,
    yaml_node_t *deme,
    struct toplevel_defaults *toplevel_defaults,
    struct demelevel_defaults *demelevel_defaults)
{
    int ret = 0;
    yaml_node_t *defaults;

    init_demelevel_defaults(demelevel_defaults);

    if ((defaults = get_value(document, deme, "defaults"))) {
        if ((ret = check_allowed_demelevel_defaults(document, defaults))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_defaults_epoch(
                document, defaults, &demelevel_defaults->epoch))) {
            goto err0;
        }
    }

    // Merge toplevel defaults.

    if (isnan(demelevel_defaults->epoch.end_time)) {
        demelevel_defaults->epoch.end_time = toplevel_defaults->epoch.end_time;
    }
    if (isnan(demelevel_defaults->epoch.start_size)) {
        demelevel_defaults->epoch.start_size = toplevel_defaults->epoch.start_size;
    }
    if (isnan(demelevel_defaults->epoch.end_size)) {
        demelevel_defaults->epoch.end_size = toplevel_defaults->epoch.end_size;
    }
    if (demelevel_defaults->epoch.size_function == -1) {
        demelevel_defaults->epoch.size_function =
            toplevel_defaults->epoch.size_function;
    }
    if (isnan(demelevel_defaults->epoch.selfing_rate)) {
        demelevel_defaults->epoch.selfing_rate =
            toplevel_defaults->epoch.selfing_rate;
    }
    if (isnan(demelevel_defaults->epoch.cloning_rate)) {
        demelevel_defaults->epoch.cloning_rate =
            toplevel_defaults->epoch.cloning_rate;
    }

err0:
    return ret;
}

/**
 * Parse the defaults node in the yaml @p document,
 * populating @p toplevel_defaults.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_toplevel_defaults(
    yaml_document_t *document,
    yaml_node_t *root,
    struct toplevel_defaults *toplevel_defaults)
{
    int ret = 0;
    yaml_node_t *defaults;

    init_toplevel_defaults(toplevel_defaults);

    if ((defaults = get_value(document, root, "defaults"))) {
        if ((ret = check_allowed_toplevel_defaults(document, defaults))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_defaults_epoch(
                document, defaults, &toplevel_defaults->epoch))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_defaults_deme(
                document, defaults, &toplevel_defaults->deme))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_defaults_migration(
                document, defaults, &toplevel_defaults->migration))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_defaults_pulse(
                document, defaults, &toplevel_defaults->pulse))) {
            goto err0;
        }
    }

err0:
    return ret;
}

static void
free_toplevel_defaults(struct toplevel_defaults *toplevel_defaults)
{
    if (toplevel_defaults->deme.ancestors) {
        free(toplevel_defaults->deme.ancestors);
    }
    if (toplevel_defaults->deme.proportions) {
        free(toplevel_defaults->deme.proportions);
    }
    if (toplevel_defaults->migration.demes) {
        free(toplevel_defaults->migration.demes);
    }
    if (toplevel_defaults->pulse.sources) {
        free(toplevel_defaults->pulse.sources);
    }
    if (toplevel_defaults->pulse.proportions) {
        free(toplevel_defaults->pulse.proportions);
    }
}

/**
 * Parse the demes node (@p demes) in the yaml @p document,
 * populating @p graph->demes.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_demes(
    yaml_document_t *document,
    struct demes_graph *graph,
    yaml_node_t *demes,
    struct toplevel_defaults *toplevel_defaults)
{
    demes_char_t **ancestors = NULL;
    double *proportions = NULL;
    int ret = 0;
    int i;
    int n_demes;

    if (demes->type != YAML_SEQUENCE_NODE) {
        errmsg("line %ld: demes must be a list of one or more deme objects\n",
                demes->start_mark.line);
        ret = DEMES_ERR_TYPE;
        goto err0;
    }

    n_demes = demes->data.sequence.items.top - demes->data.sequence.items.start;

    if (n_demes == 0) {
        errmsg("line %ld: demes must be a list of one or more deme objects\n",
                demes->start_mark.line);
        ret = DEMES_ERR_MISSING_REQUIRED;
        goto err0;
    }

    graph->demes = calloc(n_demes, sizeof *graph->demes);
    if (graph->demes == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for (i=0; i<n_demes; i++) {
        demes_char_t *name, *description;
        double start_time;
        size_t n_ancestors, n_proportions;
        yaml_node_t *deme, *epochs;
        struct demelevel_defaults demelevel_defaults;

        deme = yaml_document_get_node(document, demes->data.sequence.items.start[i]);
        assert(deme != NULL);
        if (deme->type != YAML_MAPPING_NODE) {
            errmsg("line %ld: demes[%d]: expected a deme object\n",
                    deme->start_mark.line, i);
            ret = DEMES_ERR_TYPE;
            goto err0;
        }
        if ((ret = check_allowed_demeslevel(document, deme))) {
            goto err0;
        }

        if ((ret = get_string(document, deme, "name", &name))) {
            goto err0;
        }
        if (name == NULL) {
            errmsg("line %ld: demes[%d]: required field 'name' not found\n",
                    deme->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err0;
        }

        description = toplevel_defaults->deme.description;
        start_time = toplevel_defaults->deme.start_time;
        ancestors = toplevel_defaults->deme.ancestors;
        n_ancestors = toplevel_defaults->deme.n_ancestors;
        proportions = toplevel_defaults->deme.proportions;
        n_proportions = toplevel_defaults->deme.n_proportions;

        if (get_value(document, deme, "description")) {
            if ((ret = get_string(document, deme, "description", &description))) {
                goto err0;
            }
        }

        if (get_value(document, deme, "start_time")) {
            if ((ret = get_number(document, deme, "start_time", &start_time))) {
                goto err0;
            }
        }

        if (get_value(document, deme, "ancestors")) {
            if ((ret = get_string_list(document, deme, "ancestors", &ancestors,
                            &n_ancestors))) {
                goto err0;
            }
        }

        if (get_value(document, deme, "proportions")) {
            if ((ret = get_number_list(document, deme, "proportions", &proportions,
                            &n_proportions))) {
                goto err0;
            }
        }

        if (isnan(start_time)) {
            if (n_ancestors > 0) {
                if (n_ancestors == 1) {
                    struct demes_deme *anc = demes_graph_get_deme(graph, ancestors[0]);
                    if (anc == NULL) {
                        errmsg("deme %s: ancestor deme '%s' not found. "
                                "Note: ancestor demes must be specified "
                                "before their children.\n",
                                name, ancestors[0]);
                        ret = DEMES_ERR_VALUE;
                        goto err0;
                    }
                    start_time = anc->epochs[anc->n_epochs - 1].end_time;
                } else {
                    errmsg("line %ld: deme %s: start_time is required for "
                            "demes with multiple ancestors\n",
                            deme->start_mark.line, name);
                    ret = DEMES_ERR_MISSING_REQUIRED;
                    goto err0;
                }
            } else {
                start_time = INFINITY;
            }
        }

        if ((ret = demes_graph_add_deme(
                        graph, name, description, start_time, ancestors,
                        n_ancestors, proportions, n_proportions))) {
            goto err0;
        }

        if ((ret = demes_graph_parse_demelevel_defaults(
                        document, deme, toplevel_defaults,
                        &demelevel_defaults))) {
            goto err0;
        }

        epochs = get_value(document, deme, "epochs");

        if ((ret = demes_graph_parse_epochs(
                        document, graph, &graph->demes[graph->n_demes - 1],
                        &demelevel_defaults.epoch,
                        epochs)) != 0) {
            goto err0;
        }

        if (ancestors && ancestors != toplevel_defaults->deme.ancestors) {
            free(ancestors);
            ancestors = NULL;
        }
        if (proportions && proportions != toplevel_defaults->deme.proportions) {
            free(proportions);
            proportions = NULL;
        }
    }

err0:
    if (ancestors && ancestors != toplevel_defaults->deme.ancestors) {
        free(ancestors);
    }
    if (proportions && proportions != toplevel_defaults->deme.proportions) {
        free(proportions);
    }
    return ret;
}

/**
 * Insert @p num into @p _list, maintaining value-descending sort order.
 *
 * If @p num already exists in @p _list, no action is taken.
 *
 * @param[in] num The number to be added to the list.
 * @param[in,out] _list A reference to the list of sorted numbers. This must
 *                      initially be a reference to NULL (indicating an empty
 *                      list).
 * @param[in,out] _length A reference to the list length. This must initially
 *                        be a reference to 0 (indicating an empty list).
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
insert_sorted_unique(double num, double **_list, size_t *_length)
{
    int i, j;
    int ret = 0;
    double *list = *_list;
    size_t length = *_length;
    size_t size;

    // TODO: binary search.
    for (i=0; i<length && num<list[i]; i++)
        ;

    if (i < length && num == list[i]) {
        // not unique
        return 0;
    }

    length++;
    size = length * sizeof *list;
    if ((list = realloc(list, size)) == NULL) {
        perror("realloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for (j=length - 1; j>i; j--) {
        list[j] = list[j-1];
    }
    list[i] = num;

    *_list = list;
    *_length = length;

err0:
    return ret;
}

/**
 * Build a list of migration matrices and their corresponding end times.
 *
 * There will always be at least one migration matrix, even when the graph
 * does not define any migrations. The migration matrices are not gauranteed
 * to be unique.
 *
 * The caller is responsible for freeing both @p _mm_list and @p _end_times.
 * If an error occurs, both @p _mm_list and @p _end_times will be references
 * to NULL.
 *
 * @param[in] graph The graph containing the migration matrices.
 * @param[out] _mm_list A reference to the list of migration matrices.
 *                      The migration matrices are sorted from oldest to
 *                      youngest. For m matrices and n demes, the list is
 *                      a pointer to a flat block of memory with m*n*n entries.
 *                      Entry j*n*n + k*n + l is the migration rate from
 *                      source l to dest k in matrix j.
 * @param[out] _end_times A reference to the list of end times corresponding to
 *                        each migration matrix. The end times are sorted from
 *                        oldest to youngest, and the start time for the first
 *                        migration matrix is infinity.
 * @param[out] _n_matrices A reference to the number of migration matrices.
 *                         This is also the number of end times.
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
migration_matrices(
    struct demes_graph *graph,
    double **_mm_list, double **_end_times, size_t *_n_matrices)
{
    int i, j;
    int ret = 0;
    size_t n_end_times = 0;
    size_t nmemb;
    double start_time, end_time;
    double *end_times = NULL;
    double *mm_list;

    *_mm_list = NULL;
    *_end_times = NULL;
    *_n_matrices = 0;

    /*
     * Construct the list of end times, which are the unique migration
     * start and/or end times, sorted from oldest to youngest.
     */
    for(i=0; i<graph->n_migrations; i++) {
        struct demes_migration *migration = graph->migrations + i;
        if ((ret = insert_sorted_unique(migration->end_time,
                        &end_times, &n_end_times))) {
            goto err0;
        }
        if (!isinf(migration->start_time) &&
                (ret = insert_sorted_unique(migration->start_time,
                        &end_times, &n_end_times))) {
            goto err0;
        }
    }

    // Ensure there's always 1 matrix, even when there are no migrations.
    if (n_end_times == 0) {
        int k;
        double youngest_end_time = INFINITY;
        for (k=0; k<graph->n_demes; k++) {
            struct demes_deme *deme = graph->demes + k;
            double end_time = deme->epochs[deme->n_epochs - 1].end_time;
            if (end_time < youngest_end_time) {
                youngest_end_time = end_time;
            }
        }
        insert_sorted_unique(youngest_end_time, &end_times, &n_end_times);
    }

    nmemb = n_end_times * graph->n_demes * graph->n_demes;
    if ((mm_list = calloc(nmemb, sizeof *mm_list)) == NULL) {
        perror("calloc");
        ret = DEMES_ERR_MEMORY;
        goto err0;
    }

    for(i=0; i<graph->n_migrations; i++) {
        struct demes_migration *migration = graph->migrations + i;

        start_time = INFINITY;
        for (j=0; j<n_end_times; j++) {
            end_time = end_times[j];
            if (start_time <= migration->end_time) {
                break;
            }
            if (end_time < migration->start_time) {
                int source_id = get_deme_id(graph, migration->source->name);
                int dest_id = get_deme_id(graph, migration->dest->name);
                size_t n = graph->n_demes;
                int x = j * n * n + dest_id * n + source_id;

                assert(source_id >= 0);
                assert(dest_id >= 0);

                if (mm_list[x] > 0) {
                    errmsg("multiple migrations defined for source deme %s "
                            "and dest deme %s for the time interval [%lf, %lf]\n",
                            migration->source->name, migration->dest->name,
                            start_time, end_time);
                    ret = DEMES_ERR_VALUE;
                    goto err1;
                }
                mm_list[x] = migration->rate;
            }
            start_time = end_time;
        }
    }

    *_mm_list = mm_list;
    *_end_times = end_times;
    *_n_matrices = n_end_times;

err1:
    if (ret && mm_list) {
        free(mm_list);
    }
err0:
    if (ret && end_times) {
        free(end_times);
    }
    return ret;
}

/**
 * Check that the sum of migration ingress rates don't exceed 1 for any deme,
 * and that there are no overlapping migrations.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
validate_migrations(struct demes_graph *graph)
{
    int i, j, k;
    int ret = 0;
    double start_time, end_time;
    double *mm_list= NULL, *end_times = NULL;
    size_t n_matrices;
    size_t n = graph->n_demes;

    if ((ret = migration_matrices(graph, &mm_list, &end_times, &n_matrices))) {
        goto err0;
    }

    start_time = INFINITY;
    for (i=0; i<n_matrices; i++) {
        end_time = end_times[i];
        for (j=0; j<n; j++) {
            double sum = 0.0;
            for (k=0; k<n; k++) {
                int x = i * n * n + j * n + k;
                sum += mm_list[x];
            }
            if (sum > 1.0) {
                errmsg("sum of migration rates into deme %s is greater than 1 "
                        "for the time interval (%lf, %lf]\n",
                        graph->demes[j].name, start_time, end_time);
                ret = DEMES_ERR_VALUE;
                goto err1;
            }
        }
        start_time = end_time;
    }

err1:
    free(mm_list);
    free(end_times);
err0:
    return ret;
}

/**
 * Parse the migrations node (@p migrations) in the yaml @p document,
 * populating @p graph->migrations.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_migrations(
    yaml_document_t *document,
    struct demes_graph *graph,
    yaml_node_t *migrations,
    struct toplevel_defaults *toplevel_defaults)
{
    demes_char_t **demes = NULL;
    int n_migrations;
    int ret = 0;
    int i;

    if (migrations->type != YAML_SEQUENCE_NODE) {
        errmsg("line %ld: migrations must be a list of migration objects\n",
                migrations->start_mark.line);
        ret = DEMES_ERR_TYPE;
        goto err0;
    }
    n_migrations = migrations->data.sequence.items.top
                        - migrations->data.sequence.items.start;

    for (i=0; i<n_migrations; i++) {
        yaml_node_t *migration;
        demes_char_t *source, *dest;
        size_t n_demes;
        double start_time, end_time, rate;

        migration = yaml_document_get_node(
                document, migrations->data.sequence.items.start[i]);
        assert(migration != NULL);

        if (migration->type != YAML_MAPPING_NODE) {
            errmsg("line %ld: migration[%d]: expected a migration object\n",
                    migration->start_mark.line, i);
            ret = DEMES_ERR_TYPE;
            goto err1;
        }

        if ((ret = check_allowed_migration(document, migration))) {
            goto err1;
        }

        demes = toplevel_defaults->migration.demes;
        n_demes = toplevel_defaults->migration.n_demes;
        source = toplevel_defaults->migration.source;
        dest = toplevel_defaults->migration.dest;
        start_time = toplevel_defaults->migration.start_time;
        end_time = toplevel_defaults->migration.end_time;
        rate = toplevel_defaults->migration.rate;

        if (get_value(document, migration, "demes")) {
            if ((ret = get_string_list(document, migration, "demes", &demes,
                            &n_demes))) {
                goto err1;
            }
        }
        if (get_value(document, migration, "source")) {
            if ((ret = get_string(document, migration, "source", &source))) {
                goto err1;
            }
        }
        if (get_value(document, migration, "dest")) {
            if ((ret = get_string(document, migration, "dest", &dest))) {
                goto err1;
            }
        }
        if (get_value(document, migration, "start_time")) {
            if ((ret = get_number(document, migration, "start_time",
                            &start_time))) {
                goto err1;
            }
        }
        if (get_value(document, migration, "end_time")) {
            if ((ret = get_number(document, migration, "end_time", &end_time))) {
                goto err1;
            }
        }
        if (get_value(document, migration, "rate")) {
            if ((ret = get_number(document, migration, "rate", &rate))) {
                goto err1;
            }
        }
        if (isnan(rate)) {
            errmsg("line %ld: migration[%d]: required field 'rate' not found\n",
                    migration->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err1;
        }

        if (!((source && dest && demes == NULL) ||
                    (source == NULL && dest == NULL && demes))) {
            errmsg("line %ld: migration[%d]: must specify both source and dest "
                    "demes, or alternately specify a list of deme names\n",
                    migration->start_mark.line, i);
            ret = DEMES_ERR_VALUE;
            goto err1;
        }

        if (demes) {
            if ((ret = demes_graph_add_symmetric_migration(
                            graph, demes, n_demes, start_time, end_time, rate))) {
                goto err1;
            }

            if (demes != toplevel_defaults->migration.demes) {
                free(demes);
                demes = NULL;
            }
        } else {
            if ((ret = demes_graph_add_migration(graph, source, dest,
                            start_time, end_time, rate))) {
                goto err1;
            }
        }
    }

    if ((ret = validate_migrations(graph))) {
        goto err1;
    }

err1:
    if (demes && demes != toplevel_defaults->migration.demes) {
        free(demes);
    }
err0:
    return ret;
}

/**
 * Parse the pulses node (@p pulses) in the yaml @p document,
 * populating @p toplevel_defaults->pulse.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_pulses(
    yaml_document_t *document,
    struct demes_graph *graph,
    yaml_node_t *pulses,
    struct toplevel_defaults *toplevel_defaults)
{
    demes_char_t **sources = NULL;
    double *proportions = NULL;
    int n_pulses;
    int ret = 0;
    int i;

    if (pulses->type != YAML_SEQUENCE_NODE) {
        errmsg("line %ld: pulses must be a list of pulse objects\n",
                pulses->start_mark.line);
        ret = DEMES_ERR_TYPE;
        goto err0;
    }
    n_pulses = pulses->data.sequence.items.top
                        - pulses->data.sequence.items.start;

    for (i=0; i<n_pulses; i++) {
        yaml_node_t *pulse;
        demes_char_t *dest;
        double time;
        size_t n_sources, n_proportions;

        pulse = yaml_document_get_node(
                document, pulses->data.sequence.items.start[i]);
        assert(pulse != NULL);

        if (pulse->type != YAML_MAPPING_NODE) {
            errmsg("line %ld: pulse[%d]: expected a pulse object\n",
                    pulse->start_mark.line, i);
            ret = DEMES_ERR_TYPE;
            goto err0;
        }

        if ((ret = check_allowed_pulse(document, pulse))) {
            goto err0;
        }

        sources = toplevel_defaults->pulse.sources;
        n_sources = toplevel_defaults->pulse.n_sources;
        dest = toplevel_defaults->pulse.dest;
        time = toplevel_defaults->pulse.time;
        proportions = toplevel_defaults->pulse.proportions;
        n_proportions = toplevel_defaults->pulse.n_proportions;

        if (get_value(document, pulse, "sources")) {
            if ((ret = get_string_list(document, pulse, "sources", &sources,
                            &n_sources))) {
                goto err0;
            }
            if (sources == NULL) {
                errmsg("line %ld: pulse[%d]: 'sources' must be a non-empty "
                        "list of deme names\n",
                        pulse->start_mark.line, i);
                ret = DEMES_ERR_MISSING_REQUIRED;
                goto err0;
            }
        }
        if (sources == NULL) {
            errmsg("line %ld: pulse[%d]: required field 'sources' not found\n",
                    pulse->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err0;
        }

        if (get_value(document, pulse, "dest")) {
            if ((ret = get_string(document, pulse, "dest", &dest))) {
                goto err0;
            }
        }
        if (dest == NULL) {
            errmsg("line %ld: pulse[%d]: required field 'dest' not found\n",
                    pulse->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err0;
        }

        if (get_value(document, pulse, "time")) {
            if ((ret = get_number(document, pulse, "time", &time))) {
                goto err0;
            }
        }
        if (isnan(time)) {
            errmsg("line %ld: pulse[%d]: required field 'time' not found\n",
                    pulse->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err0;
        }

        if (get_value(document, pulse, "proportions")) {
            if ((ret = get_number_list(document, pulse, "proportions",
                            &proportions, &n_proportions))) {
                goto err0;
            }
        }
        if (proportions == NULL) {
            errmsg("line %ld: pulse[%d]: required field 'proportions' not found\n",
                    pulse->start_mark.line, i);
            ret = DEMES_ERR_MISSING_REQUIRED;
            goto err0;
        }

        if ((ret = demes_graph_add_pulse(graph, sources, n_sources, dest, time,
                        proportions, n_proportions))) {
            goto err0;
        }

        if (sources && sources != toplevel_defaults->pulse.sources) {
            free(sources);
            sources = NULL;
        }
        if (proportions && proportions != toplevel_defaults->pulse.proportions) {
            free(proportions);
            proportions = NULL;
        }
    }

err0:
    if (sources && sources != toplevel_defaults->pulse.sources) {
        free(sources);
    }
    if (proportions && proportions != toplevel_defaults->pulse.proportions) {
        free(proportions);
    }
    return ret;
}

/**
 * Parse the toplevel metadata node (@p metadata) in the yaml @p document.
 * TODO: store the metadata, so we can at least round-trip it.
 *
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
static int
demes_graph_parse_toplevel_metadata(
    yaml_document_t *document,
    struct demes_graph *graph,
    yaml_node_t *metadata)
{
    int ret = 0;

    if (metadata->type != YAML_MAPPING_NODE) {
        errmsg("line %ld: metadata must be a mapping of key/value pairs\n",
                metadata->start_mark.line);
        ret = DEMES_ERR_TYPE;
        goto err0;
    }

err0:
    return ret;
}

int
demes_graph_parse(yaml_parser_t * parser, struct demes_graph **_graph)
{
    yaml_document_t document;
    yaml_node_t *root, *demes, *migrations, *pulses, *metadata;
    struct demes_graph *graph;
    struct toplevel_defaults toplevel_defaults;
    int ret = 0;
    demes_char_t *description, *time_units;
    double generation_time;
    demes_char_t **doi = NULL;
    size_t n_dois;
    locale_t locale_C, locale_saved;

    *_graph = NULL;

    /*
     * Set the C locale so that strtod() converts numbers appropriately
     * when the calling application is using another locale.
     */
    if ((locale_C = newlocale(LC_ALL, "C", (locale_t)0)) == (locale_t)0) {
        perror("newlocale");
        ret = DEMES_ERR_LOCALE;
        goto err0;
    }
    if ((locale_saved = uselocale(locale_C)) == (locale_t)0) {
        perror("uselocale");
        ret = DEMES_ERR_LOCALE;
        goto err1;
    }

    if (yaml_parser_load(parser, &document) == 0) {
        ret = DEMES_ERR_YAML;
        goto err2;
    }

    if ((root = yaml_document_get_root_node(&document)) == NULL) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    if (root->type != YAML_MAPPING_NODE) {
        ret = DEMES_ERR_TYPE;
        goto err3;
    }

    ret = check_allowed_toplevel(&document, root);
    if (ret != 0) {
        goto err3;
    }

    if ((ret = get_string(&document, root, "description", &description))) {
        goto err3;
    }
    if ((ret = get_string(&document, root, "time_units", &time_units))) {
        goto err3;
    }
    if ((ret = get_number(&document, root, "generation_time", &generation_time))) {
        goto err3;
    }

    if ((ret = get_string_list(&document, root, "doi", &doi, &n_dois))) {
        goto err3;
    }

    if ((ret = demes_graph_parse_toplevel_defaults(&document, root,
                    &toplevel_defaults))) {
        goto err4;
    }

    if ((ret = demes_graph_init(&graph, description, time_units,
                    generation_time, doi, n_dois))) {
        goto err4;
    }

    if ((demes = get_value(&document, root, "demes")) == NULL) {
        errmsg("missing demes list\n");
        ret = DEMES_ERR_MISSING_REQUIRED;
        goto err5;
    }

    if ((ret = demes_graph_parse_demes(
                    &document, graph, demes, &toplevel_defaults))) {
        goto err5;
    }

    if ((migrations = get_value(&document, root, "migrations"))) {
        if ((ret = demes_graph_parse_migrations(
                        &document, graph, migrations, &toplevel_defaults))) {
            goto err5;
        }
    }
    if ((pulses = get_value(&document, root, "pulses"))) {
        if ((ret = demes_graph_parse_pulses(
                        &document, graph, pulses, &toplevel_defaults))) {
            goto err5;
        }
    }

    if ((metadata = get_value(&document, root, "metadata"))) {
        if ((ret = demes_graph_parse_toplevel_metadata(&document, graph, metadata))) {
            goto err5;
        }
    }

    *_graph = graph;

err5:
    if (ret != 0) {
        demes_graph_free(graph);
    }
err4:
    free_toplevel_defaults(&toplevel_defaults);
err3:
    if (doi) {
        free(doi);
    }
    yaml_document_delete(&document);
err2:
    if (uselocale(locale_saved) == (locale_t)0) {
        perror("uselocale");
        if (ret == 0) {
            ret = DEMES_ERR_LOCALE;
        }
    }
err1:
    freelocale(locale_C);
err0:
    return ret;
}

int
demes_graph_load(const char *filename, struct demes_graph **graph)
{
    yaml_parser_t parser;
    FILE *fp;
    int ret;

    *graph = NULL;

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        errmsg("%s: %s\n", filename, strerror(errno));
        ret = DEMES_ERR_IO;
        goto err0;
    }

    if (yaml_parser_initialize(&parser) == 0) {
        errmsg("libyaml failed to initialise parser\n");
        ret = DEMES_ERR_YAML;
        goto err1;
    }

    yaml_parser_set_input_file(&parser, fp);
    ret = demes_graph_parse(&parser, graph);
    yaml_parser_delete(&parser);

err1:
    (void)fclose(fp);
err0:
    return ret;
}

/*
 * Functions for output.
 */

static int
append_mapping_string(
    yaml_document_t *document, int node, char *key, demes_char_t *value, yaml_scalar_style_t style)
{
    int k, v;
    int ret = 0;

    if ((k = yaml_document_add_scalar(document, NULL, (demes_char_t *)key, -1, 0)) == 0) {
        errmsg("libyaml failed to add scalar to the document (key)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

    if ((v = yaml_document_add_scalar(document, NULL, value, -1, style)) == 0) {
        errmsg("libyaml failed to add scalar to the document (string value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_mapping_pair(document, node, k, v) == 0) {
        errmsg("libyaml failed to add mapping pair to the document (string)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

err0:
    return ret;
}

static void
double_to_string(char *buf, size_t buflen, double num)
{
    assert(!isnan(num));

    if (isinf(num)) {
        assert(num > 0); // negative infinity is not valid in any field
        if (num > 0) {
            snprintf(buf, buflen, ".inf");
        } /* else {
            snprintf(buf, buflen, "-.inf");
        } */
    } else {
        /* Output sufficiently many digits to safely round-trip the
         * number. I.e., if we read the emitted YAML back in, we should
         * obtain an identical value. The YAML spec doesn't specify a
         * precision, so the value may differ in other parsers,
         * but at least we'll get consistency for testing.
         */
        snprintf(buf, buflen, "%.*lg", DBL_DECIMAL_DIG - 1, num);
    }
}

static int
append_mapping_number(
    yaml_document_t *document, int node, char *key, double value)
{
    int k, v;
    int ret = 0;
    const size_t buflen = 64;
    char number[buflen];

    double_to_string(number, buflen, value);

    if ((k = yaml_document_add_scalar(document, NULL, (demes_char_t *)key, -1, 0)) == 0) {
        errmsg("libyaml failed to add scalar to the document (key)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

    if ((v = yaml_document_add_scalar(document, NULL, //(demes_char_t*)YAML_FLOAT_TAG,
                    (demes_char_t*)number, -1, 0)) == 0) {
        errmsg("libyaml failed to add scalar to the document (float value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_mapping_pair(document, node, k, v) == 0) {
        errmsg("libyaml failed to add mapping pair to the document (float)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

err0:
    return ret;
}

static int
append_mapping_sequence(
    yaml_document_t *document, int node, char *key,
    yaml_sequence_style_t style, int *value_node)
{
    int k, v;
    int ret = 0;

    *value_node = 0;

    if ((k = yaml_document_add_scalar(document, NULL, (demes_char_t *)key, -1, 0)) == 0) {
        errmsg("libyaml failed to add scalar to document (sequence key)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if ((v = yaml_document_add_sequence(document, NULL, style)) == 0) {
        errmsg("libyaml failed to add sequence to document (sequence value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_mapping_pair(document, node, k, v) == 0) {
        errmsg("libyaml failed to add mapping pair to the document (sequence)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

    *value_node = v;

err0:
    return ret;
}

static int
append_sequence_string(
    yaml_document_t *document, int sequence, demes_char_t *value, yaml_scalar_style_t style)
{
    int item;
    int ret = 0;

    if ((item = yaml_document_add_scalar(document, NULL, value, -1, style)) == 0) {
        errmsg("libyaml failed to add scalar to the document (sequence string value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_sequence_item(document, sequence, item) == 0) {
        errmsg("libyaml failed to add mapping pair to the document (string)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

err0:
    return ret;
}

static int
append_sequence_number(
    yaml_document_t *document, int sequence, double value)
{
    int item;
    int ret = 0;
    const size_t buflen = 64;
    char number[buflen];

    double_to_string(number, buflen, value);

    if ((item = yaml_document_add_scalar(document, NULL, //(demes_char_t*)YAML_FLOAT_TAG,
                    (demes_char_t*)number, -1, 0)) == 0) {
        errmsg("libyaml failed to add scalar to the document (sequence number value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_sequence_item(document, sequence, item) == 0) {
        errmsg("libyaml failed to add sequence item to the document (float)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

err0:
    return ret;
}

static int
append_sequence_mapping(
    yaml_document_t *document, int sequence,
    yaml_mapping_style_t style, int *item_node)
{
    int item;
    int ret = 0;

    *item_node = 0;

    if ((item = yaml_document_add_mapping(document, NULL, style)) == 0) {
        errmsg("libyaml failed to add mapping to the document (sequence value)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }
    if (yaml_document_append_sequence_item(document, sequence, item) == 0) {
        errmsg("libyaml failed to add sequence item to the document (mapping)\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

    *item_node = item;

err0:
    return ret;
}

int
demes_graph_emit(struct demes_graph *graph, yaml_emitter_t *emitter)
{
    locale_t locale_C, locale_saved;
    yaml_document_t document;
    int ret = 0;
    int root, doi_node, demes_node, pulses_node, migrations_node;
    int i, j;

    /*
     * Set the C locale so that we output valid YAML numbers, even when
     * the calling application has set a locale that uses a decimal comma
     * rather than a decimal point.
     */
    if ((locale_C = newlocale(LC_ALL, "C", (locale_t)0)) == (locale_t)0) {
        perror("newlocale");
        ret = DEMES_ERR_LOCALE;
        goto err0;
    }
    if ((locale_saved = uselocale(locale_C)) == (locale_t)0) {
        perror("uselocale");
        ret = DEMES_ERR_LOCALE;
        goto err1;
    }

    if (yaml_document_initialize(&document, NULL, NULL, NULL, 1, 1) == 0) {
        errmsg("libyaml failed to initialise document\n");
        ret = DEMES_ERR_YAML;
        goto err2;
    }

    if ((root = yaml_document_add_mapping(&document, NULL,
                    YAML_BLOCK_MAPPING_STYLE)) == 0) {
        errmsg("libyaml failed to create root node\n");
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    if (graph->description) {
        if ((ret = append_mapping_string(&document, root, "description",
                        graph->description, YAML_DOUBLE_QUOTED_SCALAR_STYLE))) {
            goto err3;
        }
    } else {
        if ((ret = append_mapping_string(&document, root, "description",
                        (demes_char_t*)"", YAML_DOUBLE_QUOTED_SCALAR_STYLE))) {
            goto err3;
        }
    }

    {
        yaml_scalar_style_t style = 0;
        if (strcmp((char*)graph->time_units, "generations")) {
            // User could put arbitrary text in here, so we must quote it to be
            // robust to round-tripping. (E.g. if it looks like a number)
            style = YAML_DOUBLE_QUOTED_SCALAR_STYLE;
        }
        if ((ret = append_mapping_string(&document, root, "time_units",
                        graph->time_units, style))) {
            goto err3;
        }
    }

    if ((ret = append_mapping_number(&document, root, "generation_time",
                    graph->generation_time))) {
        goto err3;
    }

    if ((ret = append_mapping_sequence(&document, root, "doi",
                    YAML_BLOCK_SEQUENCE_STYLE, &doi_node))) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }
    if (graph->doi) {
        for (i=0; i<graph->n_dois; i++) {
            if ((ret = append_sequence_string(&document, doi_node, graph->doi[i], YAML_DOUBLE_QUOTED_SCALAR_STYLE))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
        }
    }

    if ((ret = append_mapping_sequence(&document, root, "demes",
                    YAML_BLOCK_SEQUENCE_STYLE, &demes_node))) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    for (i=0; i<graph->n_demes; i++) {
        struct demes_deme *deme = graph->demes + i;
        int deme_node, epochs_node, ancestors_node, proportions_node;

        if ((ret = append_sequence_mapping(&document, demes_node,
                        YAML_BLOCK_MAPPING_STYLE, &deme_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_string(&document, deme_node, "name", deme->name, 0))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if (deme->description) {
            if ((ret = append_mapping_string(&document, deme_node,
                            "description", deme->description,
                            YAML_DOUBLE_QUOTED_SCALAR_STYLE))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
        } else {
            if ((ret = append_mapping_string(&document, deme_node,
                            "description", (demes_char_t*)"",
                            YAML_DOUBLE_QUOTED_SCALAR_STYLE))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
        }
        if ((ret = append_mapping_number(&document, deme_node, "start_time",
                        deme->start_time))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_sequence(&document, deme_node,
                        "ancestors", YAML_FLOW_SEQUENCE_STYLE,
                        &ancestors_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_sequence(&document, deme_node,
                        "proportions", YAML_FLOW_SEQUENCE_STYLE,
                        &proportions_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if (deme->ancestors) {
            for (j=0; j<deme->n_ancestors; j++) {
                if ((ret = append_sequence_string(&document, ancestors_node,
                                deme->ancestors[j]->name, 0))) {
                    ret = DEMES_ERR_YAML;
                    goto err3;
                }
                if ((ret = append_sequence_number(&document, proportions_node,
                                deme->proportions[j]))) {
                    ret = DEMES_ERR_YAML;
                    goto err3;
                }
            }
        }

        if ((ret = append_mapping_sequence(&document, deme_node,
                        "epochs", YAML_BLOCK_SEQUENCE_STYLE,
                        &epochs_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        for (j=0; j<deme->n_epochs; j++) {
            struct demes_epoch *epoch = deme->epochs + j;
            int epoch_node;

            if ((ret = append_sequence_mapping(&document, epochs_node,
                            YAML_BLOCK_MAPPING_STYLE, &epoch_node))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }

            if ((ret = append_mapping_number(&document, epoch_node, "end_time",
                            epoch->end_time))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
            if ((ret = append_mapping_number(&document, epoch_node, "start_size",
                            epoch->start_size))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
            if ((ret = append_mapping_number(&document, epoch_node, "end_size",
                            epoch->end_size))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }

            assert(epoch->size_function == DEMES_SIZE_FUNCTION_CONSTANT ||
                   epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL);
            if (epoch->size_function == DEMES_SIZE_FUNCTION_CONSTANT) {
                if ((ret = append_mapping_string(&document, epoch_node,
                                "size_function", (demes_char_t*)"constant", 0))) {
                    ret = DEMES_ERR_YAML;
                    goto err3;
                }
            } else if (epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL) {
                if ((ret = append_mapping_string(&document, epoch_node,
                                "size_function", (demes_char_t*)"exponential", 0))) {
                    ret = DEMES_ERR_YAML;
                    goto err3;
                }
            } /* else {
                errmsg("unknown size_function: %d\n", epoch->size_function);
                ret = DEMES_ERR_VALUE;
                goto err3;
            } */

            if ((ret = append_mapping_number(&document, epoch_node, "selfing_rate",
                            epoch->selfing_rate))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
            if ((ret = append_mapping_number(&document, epoch_node, "cloning_rate",
                            epoch->cloning_rate))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
        }
    }

    if ((ret = append_mapping_sequence(&document, root, "migrations",
                    YAML_BLOCK_SEQUENCE_STYLE, &migrations_node))) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    for (i=0; i<graph->n_migrations; i++) {
        struct demes_migration *migration = graph->migrations + i;
        int migration_node;

        if ((ret = append_sequence_mapping(&document, migrations_node,
                        YAML_FLOW_MAPPING_STYLE, &migration_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_string(&document, migration_node, "source",
                        migration->source->name, 0))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_string(&document, migration_node, "dest",
                        migration->dest->name, 0))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_number(&document, migration_node, "start_time",
                        migration->start_time))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_number(&document, migration_node, "end_time",
                        migration->end_time))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_number(&document, migration_node, "rate",
                        migration->rate))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
    }

    if ((ret = append_mapping_sequence(&document, root, "pulses",
                    YAML_BLOCK_SEQUENCE_STYLE, &pulses_node))) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    for (i=0; i<graph->n_pulses; i++) {
        struct demes_pulse *pulse = graph->pulses + i;
        int pulse_node, sources_node, proportions_node;

        if ((ret = append_sequence_mapping(&document, pulses_node,
                        YAML_FLOW_MAPPING_STYLE, &pulse_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_string(&document, pulse_node, "dest",
                        pulse->dest->name, 0))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_number(&document, pulse_node, "time",
                        pulse->time))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }

        if ((ret = append_mapping_sequence(&document, pulse_node,
                        "sources", YAML_FLOW_SEQUENCE_STYLE,
                        &sources_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        if ((ret = append_mapping_sequence(&document, pulse_node,
                        "proportions", YAML_FLOW_SEQUENCE_STYLE,
                        &proportions_node))) {
            ret = DEMES_ERR_YAML;
            goto err3;
        }
        for (j=0; j<pulse->n_sources; j++) {
            if ((ret = append_sequence_string(&document, sources_node,
                            pulse->sources[j]->name, 0))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
            if ((ret = append_sequence_number(&document, proportions_node,
                            pulse->proportions[j]))) {
                ret = DEMES_ERR_YAML;
                goto err3;
            }
        }
    }

    if (yaml_emitter_open(emitter) == 0) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    if (yaml_emitter_dump(emitter, &document) == 0) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

    if (yaml_emitter_close(emitter) == 0) {
        ret = DEMES_ERR_YAML;
        goto err3;
    }

err3:
    yaml_document_delete(&document);
err2:
    if (uselocale(locale_saved) == (locale_t)0) {
        perror("uselocale");
        if (ret == 0) {
            ret = DEMES_ERR_LOCALE;
        }
    }
err1:
    freelocale(locale_C);
err0:
    return ret;
}

int
demes_graph_dump(struct demes_graph *graph, FILE *fp)
{
    yaml_emitter_t emitter;
    int ret;

    if (yaml_emitter_initialize(&emitter) == 0) {
        errmsg("libyaml failed to initialise emitter\n");
        ret = DEMES_ERR_YAML;
        goto err0;
    }

    yaml_emitter_set_output_file(&emitter, fp);
    ret = demes_graph_emit(graph, &emitter);
    yaml_emitter_delete(&emitter);

err0:
    return ret;
}
