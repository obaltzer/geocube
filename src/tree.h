/**
 * @file tree.h 
 *
 * Defines functions to create main memory and external memory query trees.
 */
#ifndef __TREE_H
#define __TREE_H

#include "hilbert.h"
#include "sort.h"

/**
 * Enumeration of node types in the tree.
 */
enum tnode_type
{
    INTERMEDIATE,
    LEAF
};

/**
 * Structure representing a node in an internal memory tree. The children
 * pointer is interpreted depending on the type of the node. An instance of
 * this struct is always followed by a memory fragment that represents the
 * bounding box of the node.
 */
struct fp_im_tnode
{
    enum tnode_type type;
    size_t level;
    size_t n_children;
    void* children;
    size_t n_leaves;
    void* leaves;
};

struct fp_im_tnode* fp_im_build_tree(struct fp_context* context,
                                     size_t n, void* input, 
                                     enum tnode_type type);
void fp_im_destroy_tree(struct fp_im_tnode* tree);
void fp_print_bbox(struct fp_context* context, struct fp_im_tnode* node);
size_t fp_im_query_tree(struct fp_context* context, struct fp_im_tnode* root,
                     void* min, void* max, 
                     void (*callback)(struct fp_context*, size_t, void*));
#endif
