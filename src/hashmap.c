/*****************************************************************************
 * This file is part of uvg266 VVC encoder.
 *
 * Copyright (c) 2023, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "hashmap.h"

/**
 * \brief This function creates a node for the uvg_hashmap.
 * 
 * \param key   the key of the node to be created
 * \param value the value of the node to be created
 * \return uvg_hashmap_node a node with the given key and value
 */
uvg_hashmap_node_t* uvg_hashmap_create_node(uint32_t key, uint32_t value) {
  uvg_hashmap_node_t* new_node = (uvg_hashmap_node_t*)malloc(sizeof(uvg_hashmap_node_t));
  new_node->key = key;
  new_node->value = value;
  new_node->next = NULL;
  new_node->size  = 1;
  return new_node;
}

/**
 * \brief This function creates a new uvg_hashmap with a given bucket size.
 * 
 * \param bucket_size the size of the hashmap bucket
 * \return uvg_hashmap a new uvg_hashmap with the given bucket size
 */
uvg_hashmap_t* uvg_hashmap_create(uint32_t bucket_size)
{
  uvg_hashmap_t* new_hashmap = (uvg_hashmap_t*)malloc(sizeof(uvg_hashmap_t));
  new_hashmap->bucket_size = bucket_size;
  new_hashmap->table = (uvg_hashmap_node_t**)malloc(sizeof(uvg_hashmap_node_t*) * bucket_size);
  for (unsigned i = 0; i < bucket_size; i++) {
    new_hashmap->table[i] = NULL;
  }
  return new_hashmap;
}

/**
 * \brief This function calculates the hash index for a given 
 *        key and bucket size using the Jenkins hash function.
 *
 * \param key         the key to be hashed
 * \param bucket_size the size of the hashmap bucket
 * \return the hashed index for the given key and bucket size. 
 */
static uint32_t uvg_hashmap_hash(uint32_t key, uint32_t bucket_size)
{
  //key ^= (key >> 20) ^ (key >> 12);
  //return (key ^ (key >> 7) ^ (key >> 4) ^ 2654435769U) % bucket_size;
  return key % bucket_size;
}

/**
 * \brief This function inserts a new node into the hashmap.
 * 
 * \param map   the hashmap to insert the new node into
 * \param key   the key of the new node
 * \param value the value of the new node
 */
void uvg_hashmap_insert(uvg_hashmap_t* map, uint32_t key, uint32_t value) {
    uint32_t hash_index = uvg_hashmap_hash(key, map->bucket_size);
    uvg_hashmap_node_t* new_node = uvg_hashmap_create_node(key, value);
    new_node->next = (void*)map->table[hash_index];
    if (new_node->next != NULL) new_node->size = ((uvg_hashmap_node_t*)new_node->next)->size + 1;
    map->table[hash_index] = new_node;
}

/**
 * \brief This function searches the hashmap for the given key.
 * 
 * \param map the hashmap to search in
 * \param key the key to search for
 * \return uvg_hashmap_node the node with the given key, NULL if not found.
 */
uvg_hashmap_node_t* uvg_hashmap_search(uvg_hashmap_t* map, uint32_t key) {
  uint32_t hashIndex = uvg_hashmap_hash(key, map->bucket_size);
  return map->table[hashIndex];
}

uint32_t uvg_hashmap_search_return_first(uvg_hashmap_t* map, uint32_t key)
{
  uint32_t            hashIndex   = uvg_hashmap_hash(key, map->bucket_size);
  uvg_hashmap_node_t* temp        = map->table[hashIndex];  
  // Search key in chain and return the first match
  while (temp) {
    if (temp->key == key) {
      return temp->value;
    }
    temp = (uvg_hashmap_node_t*)temp->next;
  }
  return -1;
}

/**
 * \brief This function frees the memory of a given hashmap node.
 * 
 * \param node the node to free the memory of.
 */
void uvg_hashmap_node_free(uvg_hashmap_node_t* node)
{
  while (node) {
    uvg_hashmap_node_t* to_delete = node;
    node                        = (uvg_hashmap_node_t*)node->next;
    free(to_delete);
  }
}

/**
 * \brief This function frees the memory of a given hashmap.
 * 
 * \param map the hashmap to free the memory of.
 */
void uvg_hashmap_free(uvg_hashmap_t* map) {
  for (unsigned i = 0; i < map->bucket_size; i++) {
    uvg_hashmap_node_t* temp = map->table[i];
    uvg_hashmap_node_free(temp);
  }
  free(map->table);
  free(map);
}
