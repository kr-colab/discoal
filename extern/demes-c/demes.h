#ifndef _DEMES_DEMES_H
#define _DEMES_DEMES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <yaml.h>

/*
 * The character type (UTF-8 octet).
 * libyaml takes care of validating and converting YAML input to utf8 octets.
 */
/* typedef unsigned char yaml_char_t; */
typedef yaml_char_t demes_char_t;

/* error codes */
enum {
    DEMES_SUCCESS=0,
    DEMES_ERR_MEMORY, /* calloc, realloc, strdup, etc. */
    DEMES_ERR_IO, /* fopen */
    DEMES_ERR_LOCALE, /* newlocale, uselocale */
    DEMES_ERR_YAML, /* error from libyaml */
    DEMES_ERR_TYPE, /* type doesn't match the data model */
    DEMES_ERR_KEY, /* invalid field */
    DEMES_ERR_VALUE, /* invalid value */
    DEMES_ERR_MISSING_REQUIRED, /* required field is missing */
    DEMES_ERR_PROPORTIONS, /* invalid proportion(s) */
};

enum size_function {
    DEMES_SIZE_FUNCTION_CONSTANT,
    DEMES_SIZE_FUNCTION_EXPONENTIAL,
};

struct demes_epoch {
    double end_time;
    double start_size;
    double end_size;
    enum size_function size_function;
    double selfing_rate;
    double cloning_rate;
};

struct demes_deme {
    demes_char_t *name;
    demes_char_t *description;
    double start_time;
    struct demes_deme **ancestors;
    double *proportions;
    struct demes_epoch *epochs;
    size_t n_ancestors;
    size_t n_epochs;
};

/* asymmetric migration */
struct demes_migration {
    struct demes_deme *source;
    struct demes_deme *dest;
    double start_time;
    double end_time;
    double rate;
};

struct demes_pulse {
    struct demes_deme **sources;
    struct demes_deme *dest;
    double time;
    double *proportions;
    size_t n_sources;
};

struct demes_graph {
    demes_char_t *description;
    demes_char_t *time_units;
    double generation_time;
    demes_char_t **doi; // list of strings
    struct demes_deme *demes;
    struct demes_migration *migrations;
    struct demes_pulse *pulses;
    size_t n_dois;
    size_t n_demes;
    size_t n_migrations;
    size_t n_pulses;
};

/**
 * Load and resolve a demes @p graph from the given libyaml @p parser.
 *
 * If successful, the address of the newly allocated graph will be returned
 * to the caller in the @p graph parameter. It's the responsibility of the
 * caller to free the memory used by the graph with demes_graph_free().
 *
 * @code{.c}
 *
 *      struct demes_graph *graph;
 *      yaml_parser_t parser;
 *      int ret;
 *
 *      if (yaml_parser_initialize(&parser) == 0) {
 *          fprintf(stderr, "failed to init yaml parser\n");
 *          exit(EXIT_FAILURE);
 *      }
 *
 *      ... // set the input for the yaml parser (see yaml.h)
 *
 *      ret = demes_graph_parse(&parser, &graph);
 *      yaml_parser_delete(&parser);
 *      if (ret) {
 *          fprintf(stderr, "failed to resolve the graph\n");
 *          exit(EXIT_FAILURE);
 *      }
 *
 *     ... // do something with the graph
 *
 *     demes_graph_free(graph);
 *  
 * @endcode
 *
 * @param[in] parser The libyaml parser from which the graph will be loaded.
 * @param[out] graph The address of a struct demes_graph pointer that will
 *                   will point to the fully-resolved graph upon successful
 *                   return from the function.
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
int demes_graph_parse(yaml_parser_t *parser, struct demes_graph **graph);

/**
 * Load and resolve a demes graph from a YAML file.
 *
 * If successful, the address of the newly allocated graph will be returned
 * to the caller in the @p graph parameter. It's the responsibility of the
 * caller to free the memory used by the graph with demes_graph_free().
 *
 * @code{.c}
 *
 *     struct demes_graph *graph;
 *     int ret;
 *
 *     if ((ret = demes_graph_parse("file.yaml", &graph))) {
 *         fprintf(stderr, "failed to load file.yaml\n");
 *         exit(EXIT_FAILURE);
 *     }
 *
 *     ... // do something with the graph
 *
 *     demes_graph_free(graph);
 *
 * @endcode
 * 
 * @param[in] filename The path to the YAML file to be loaded.
 * @param[out] graph The address of a struct demes_graph pointer that will
 *                   will point to the fully-resolved graph upon successful
 *                   return from the function.
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
int demes_graph_load(const char *filename, struct demes_graph **graph);

/**
 * Write the fully-qualified @p graph to the libyaml @p emitter.
 *
 * @param[in] graph The graph.
 * @param[in] emitter The libyaml emitter.
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
int demes_graph_emit(struct demes_graph *graph, yaml_emitter_t *emitter);

/**
 * Write the fully-qualified @p graph as YAML to the file pointer @p fp.
 *
 * @param[in] graph The graph.
 * @param[in] fp The file pointer to write to.
 * @return 0 upon success, or one of the DEMES_ERR_* error codes upon error.
 */
int demes_graph_dump(struct demes_graph *graph, FILE *fp);

/**
 * Free the memory used by the @p graph.
 *
 * @param[in] graph The graph to be freed.
 */
void demes_graph_free(struct demes_graph *graph);

#ifdef __cplusplus
}
#endif

#endif /* _DEMES_DEMES_H */
