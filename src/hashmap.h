#pragma once

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// The ratio of the hashmap bucket size to the maximum number of elements
#define UVG_HASHMAP_RATIO 12.0
// Use Hashmap for 4x4 blocks
#define UVG_HASHMAP_BLOCKSIZE 8

typedef struct uvg_hashmap_node {
    void*    next;  
    uint32_t key;
    uint32_t value;
    uint32_t size;
} uvg_hashmap_node_t;

typedef struct uvg_hashmap {
  uint32_t bucket_size;
  uvg_hashmap_node_t** table;
} uvg_hashmap_t;

uvg_hashmap_node_t* uvg_hashmap_create_node(uint32_t key, uint32_t value);

uvg_hashmap_t* uvg_hashmap_create(uint32_t bucket_size);

//uint32_t uvg_hashmap_hash(uint32_t key, uint32_t bucket_size);

void uvg_hashmap_insert(uvg_hashmap_t* map, uint32_t key, uint32_t value);

uvg_hashmap_node_t* uvg_hashmap_search(uvg_hashmap_t* map, uint32_t key);

uint32_t uvg_hashmap_search_return_first(uvg_hashmap_t* map, uint32_t key);

void uvg_hashmap_node_free(uvg_hashmap_node_t* node);

void uvg_hashmap_free(uvg_hashmap_t* map);
