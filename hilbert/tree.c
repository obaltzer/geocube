#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "hilbert.h"
#include "sort.h"
#include "tree.h"

void fp_print_bbox(struct fp_context* context, struct fp_im_tnode* node)
{
    size_t i;
    /* compute the bbox array addresses */
    void* bbox_base = ((void*)node) + sizeof(struct fp_im_tnode);
    fpz_t* minz = (fpz_t*)(bbox_base + context->minz);
    fpf_t* minf = (fpf_t*)(bbox_base + context->minf);
    fpz_t* maxz = (fpz_t*)(bbox_base + context->maxz);
    fpf_t* maxf = (fpf_t*)(bbox_base + context->maxf);
     
    printf("Min: (");
    for(i = 0; i < context->env->dimz; i++)
        if(i)
            printf(", %llu", minz[i]);
        else
            printf("%llu", minz[i]);
    for(i = 0; i < context->env->dimf; i++)
        if(i + context->env->dimz)
            printf(", %f", minf[i]);
        else
            printf("%f", minf[i]);
    printf(")\tMax: (");
    for(i = 0; i < context->env->dimz; i++)
        if(i)
            printf(", %llu", maxz[i]);
        else
            printf("%llu", maxz[i]);
    for(i = 0; i < context->env->dimf; i++)
        if(i + context->env->dimz)
            printf(", %f", maxf[i]);
        else
            printf("%f", maxf[i]);
    printf(")");
}
        
void fp_set_bbox(struct fp_context* context, struct fp_im_tnode* node)
{
    size_t i, j;
    /* compute the bbox array addresses */
    void* bbox_base = ((void*)node) + sizeof(struct fp_im_tnode);
    fpz_t* minz = (fpz_t*)(bbox_base + context->minz);
    fpf_t* minf = (fpf_t*)(bbox_base + context->minf);
    fpz_t* maxz = (fpz_t*)(bbox_base + context->maxz);
    fpf_t* maxf = (fpf_t*)(bbox_base + context->maxf);
    
    if(node->type == LEAF)
    {   
        /* create bounding box for leave node */
        void* record;
        fpz_t* coordsz;
        fpf_t* coordsf;
        
        for(i = 0; i < node->n_children; i++)
        {
            /* iterate through all records underneath the current node */
            record = node->children + (context->sizeof_record * i);

            if(i == 0)
            {
                /* the first record is initialization so just copy */
                memcpy(minz, record + context->offset_coordsz,
                       context->sizeof_record - context->offset_coordsz);
                memcpy(maxz, record + context->offset_coordsz, 
                       context->sizeof_record - context->offset_coordsz);
            }
            else
            {
                coordsz = (fpz_t*)(record + context->offset_coordsz);
                coordsf = (fpf_t*)(record + context->offset_coordsf);
                /* all other records compare with current bbox */
                for(j = 0; j < context->env->dimz; j++)
                {
                    
                    if(minz[j] > coordsz[j])
                        minz[j] = coordsz[j];
                    if(maxz[j] < coordsz[j])
                        maxz[j] = coordsz[j];
                }
                for(j = 0; j < context->env->dimf; j++)
                {
                    if(minf[j] > coordsf[j])
                        minf[j] = coordsf[j];
                    if(maxf[j] < coordsf[j])
                        maxf[j] = coordsf[j];
                } 
            }
        }
    }
    else if(node->type == INTERMEDIATE)
    {
        /* the node is an intermediate node, so we generate the bounding
         * box out of the bounding boxes of the children */
        void* c_bbox;
        fpz_t* c_minz;
        fpf_t* c_minf;
        fpz_t* c_maxz;
        fpf_t* c_maxf;
        size_t node_size = sizeof(struct fp_im_tnode) + context->bbox_size;
        
        for(i = 0; i < node->n_children; i++)
        {
            c_bbox = node->children + (node_size * i) 
                        + sizeof(struct fp_im_tnode);
            c_minz = (fpz_t*)(c_bbox + context->minz);
            c_minf = (fpf_t*)(c_bbox + context->minf);
            c_maxz = (fpz_t*)(c_bbox + context->maxz);
            c_maxf = (fpf_t*)(c_bbox + context->maxf);
            
            if(i == 0)
                /* the first child is only for initialization, just copy */
                memcpy(bbox_base, c_bbox, context->bbox_size);
            else
            {
                for(j = 0; j < context->env->dimz; j++)
                {
                    if(maxz[j] < c_maxz[j])
                        maxz[j] = c_maxz[j];
                    if(minz[j] > c_minz[j])
                        minz[j] = c_minz[j];
                }
                for(j = 0; j < context->env->dimf; j++)
                {
                    if(minf[j] > c_minf[j])
                        minf[j] = c_minf[j];
                    if(maxf[j] < c_maxf[j])
                        maxf[j] = c_maxf[j];
                } 
            }
        }
    }
}        

/**
 * Builds a new internal memory tree based on the previously sorted record
 * buffer.
 *
 * @param[in] context reference to the context structure
 * @param[in] n the number of records in the buffer
 * @param[in] input the Hilbert-sorted input buffer
 * @param[in] type specifies the type of nodes that should be created
 * @return a pointer to the root node of the tree
 */
struct fp_im_tnode* fp_im_build_tree(struct fp_context* context,
                                     size_t n, void* input, 
                                     enum tnode_type type)
{
    size_t i;
    struct fp_im_tnode* node;
    struct fp_im_tnode* last_child;
    /* number of nodes necessary*/
    size_t n_nodes = n % context->fanout ? n / context->fanout + 1
                        : n / context->fanout;
    /* compute the size of a node including its bounding box */
    size_t node_size = sizeof(struct fp_im_tnode) + context->bbox_size;
    
    /* allocate memory for the nodes and bounding boxes */
    void* nodes = malloc(node_size * n_nodes);
   
    /* XXX PROFILING */
    context->build_tree_calls++;
    context->n_tree_nodes += n_nodes;

    /* for each node group the appropriate records together */
    for(i = 0; i < n_nodes; i++)
    {
        node = (struct fp_im_tnode*)(nodes + node_size * i);
            
        node->type = type;
        /* set the number of children of the node */
        node->n_children = i < n / context->fanout
                        ? context->fanout 
                        : n % context->fanout;
        /* set the start pointer for the children of this node -- use
         * records size as offset for leave nodes and node size for
         * intermediate nodes. */
        if(type == LEAF)
        {
            node->children = input + (i * context->fanout) 
                                            * context->sizeof_record;
            /* compute the number of leaves under the node */
            node->n_leaves = node->n_children;
            node->leaves = node->children;
            node->level = 1;
            assert(node->children != NULL);
        }
        else
        {
            node->children = input + (i * context->fanout * node_size);
            /* compute the number of leaves under the node */
            if(i < n_nodes - 1)
            {
                /* this is for a node with full fanout: all its children
                 * have full fanout, hence the number of the left most
                 * child times the fanout of the fanout of the current 
                 * node */
                node->n_leaves = ((struct fp_im_tnode*)input)->n_leaves
                                    * node->n_children;
            }
            else
            {
                /* this is for a node with partial fanout */
             
                /* this node's last child is the last in the input array */
                last_child = 
                    (struct fp_im_tnode*)(input + node_size * (n - 1));

                /* the number of leaves of this node are full fanout
                 * children but one, plus the number of leaves under the
                 * last child */   
                node->n_leaves =
                    ((struct fp_im_tnode*)input)->n_leaves 
                            * (node->n_children - 1)
                            + last_child->n_leaves;
            }
            /* the leaves of this node start where the leaves of the
             * left-most child starts */
            node->leaves = 
                ((struct fp_im_tnode*)node->children)->leaves;
            node->level = 
                ((struct fp_im_tnode*)node->children)->level + 1;
        }
        /* set the boundary box of this node */
        fp_set_bbox(context, node);
    }
    /* recursively build the other levels up to the root */
    return n_nodes == 1 ? nodes 
            : fp_im_build_tree(context, n_nodes, nodes, INTERMEDIATE);
}

/**
 * Frees the memory that has been allocated for the tree and thus destroys
 * the indexing tree. This function will not free the memory that is
 * allocated to hold the actual records.
 *
 * @param[in] tree pointer to the root node of the tree.
 */
void fp_im_destroy_tree(struct fp_im_tnode* tree)
{   
    if(tree->type != LEAF)
        /* as long as we do not hit a LEAF node we want to traverse down
         * the most left branch of the tree and free each level of nodes.
         */
        fp_im_destroy_tree(tree->children);
    
    /* now delete the current level (leaf or intermediate), we do not need
     * to destroy bounding boxes since they are part of this array */
    free(tree);
    return;
}

/**
 * Tests if a record is contained in the given range.
 *
 * @param[in] context reference to the computation context
 * @param[in] min lower bound of the range
 * @param[in] max higher bound of the range
 * @param[in] r reference to the record
 * @return 1 if the record is contained in the range, otherwise 0
 */
int fp_contains_record(struct fp_context* context, void* min, void* max,
                       void* r)
{
    size_t i = 0;
    for(; i < context->env->dimz; i++)
        if(((fpz_t*)(r + context->offset_coordsz))[i] 
                < ((fpz_t*)(min + context->offset_coordsz))[i] 
                || ((fpz_t*)(r + context->offset_coordsz))[i] 
                    > ((fpz_t*)(max + context->offset_coordsz))[i])
            return 0;
    for(i = 0; i < context->env->dimf; i++)
        if(((fpf_t*)(r + context->offset_coordsf))[i] 
                < ((fpf_t*)(min + context->offset_coordsf))[i] 
                || ((fpf_t*)(r + context->offset_coordsf))[i] 
                    > ((fpf_t*)(max + context->offset_coordsf))[i])
            return 0;
    return 1;
}
  
/**
 * Checks if the bounding box of the given node is contained in the range
 * specified.
 *
 * @param[in] context reference to the processing context
 * @param[in] min lower end of the range
 * @param[in] max higher end of the range
 * @param[in] r reference to the tree node
 * @return 1 if the record's bounding box is contained in the range
 *           otherwise 0
 */
int fp_contains_node(struct fp_context* context, void* min, void* max,
                     struct fp_im_tnode* node)
{
    size_t i;
    /* compute the bbox array addresses */
    void* bbox_base = ((void*)node) + sizeof(struct fp_im_tnode);
    fpz_t* minz = (fpz_t*)(bbox_base + context->minz);
    fpf_t* minf = (fpf_t*)(bbox_base + context->minf);
    fpz_t* maxz = (fpz_t*)(bbox_base + context->maxz);
    fpf_t* maxf = (fpf_t*)(bbox_base + context->maxf);

    for(i = 0; i < context->env->dimz; i++)
        if(minz[i] < ((fpz_t*)(min + context->offset_coordsz))[i]
                || maxz[i] > ((fpz_t*)(max + context->offset_coordsz))[i])
            return 0;
    for(i = 0; i < context->env->dimf; i++)
        if(minf[i] < ((fpf_t*)(min + context->offset_coordsf))[i]
                || maxf[i] > ((fpf_t*)(max + context->offset_coordsf))[i])
            return 0;
    return 1;
}

/**
 * Queries the tree in a DFS fashion and reports sets of records that are
 * contained in the query region.
 *
 * @param[in] context reference to the environment's context
 * @param[in] root the root node of the tree to query on
 * @param[in] min the minimum corner of the query region
 * @param[in] max the maximum corner of the query region
 * @param[in] callback callback function that is called to report the set
 *            of records
 * @return overall number of records reported
 */
size_t fp_im_query_tree(struct fp_context* context, 
                        struct fp_im_tnode* root,
                        void* min, void* max,
                        void (*callback)(struct fp_context*, size_t, void*))
{
    size_t n_leaves = 0;
    size_t node_size = sizeof(struct fp_im_tnode) + context->bbox_size;
    
    if(fp_contains_node(context, min, max, root))
    {
        /* if the entire node is contained in the range, report all the
         * node's leaves */
        printf("Report node!\n");
        callback(context, root->n_leaves, root->leaves);
        n_leaves = root->n_leaves;
    }
    else
    {
        size_t i = 0;
        /* if the node is not entirely contained in the range we want to
         * traverse to its children in DFS fashion. */
        if(root->type != LEAF)
        {
            /* only call recursive of the node is an intermediate node */
            for(; i < root->n_children; i++)
                n_leaves += fp_im_query_tree(
                    context, 
                    (struct fp_im_tnode*)(root->children + i * node_size),
                    min, max, callback);
        }
        else
        {
            void* record;
            /* if all the children are leaves/records do not recurse and
             * collect information from here */
            for(; i < root->n_children; i++)
            {
                record = root->children + (i * context->sizeof_record);
                assert(record != NULL);
                if(fp_contains_record(context, min, max, record))
                {
                    printf("Report single record!\n");
                    callback(context, 1, record);
                    n_leaves++;
                }
            }
        }
    }
    return n_leaves;
}
