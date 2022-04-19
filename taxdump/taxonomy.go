// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package taxdump

import (
	"bufio"
	"errors"
	"fmt"
	"strconv"
	"strings"
	"sync"

	"github.com/shenwei356/xopen"
)

// Taxonomy holds relationship of taxon in a taxonomy.
type Taxonomy struct {
	file     string
	rootNode uint32

	Nodes      map[uint32]uint32 // child -> parent
	DelNodes   map[uint32]struct{}
	MergeNodes map[uint32]uint32 // from -> to
	Names      map[uint32]string

	taxid2rankid map[uint32]uint8 // taxid -> rank id
	ranks        []string         // rank id -> rank
	Ranks        map[string]interface{}

	hasRanks      bool
	hasDelNodes   bool
	hasMergeNodes bool
	hasNames      bool

	cacheLCA bool
	lcaCache sync.Map

	maxTaxid uint32
}

// ErrIllegalColumnIndex means column index is 0 or negative.
var ErrIllegalColumnIndex = errors.New("taxdump: illegal column index, positive integer needed")

// ErrRankNotLoaded means you should reate load Taxonomy with NewTaxonomyWithRank before calling some methods.
var ErrRankNotLoaded = errors.New("taxdump: taxonomic ranks not loaded, please call: NewTaxonomyWithRank")

// ErrNamesNotLoaded means you should call LoadNames before using taxonomy names.
var ErrNamesNotLoaded = errors.New("taxdump: taxonomy names not loaded, please call: LoadNames")

// ErrTooManyRanks means number of ranks exceed limit of 255
var ErrTooManyRanks = errors.New("taxdump: number of ranks exceed limit of 255")

// ErrUnkownRank indicate an unknown rank
var ErrUnkownRank = errors.New("taxdump: unknown rank")

// NewTaxonomyFromNCBI parses nodes relationship from nodes.dmp
// from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz .
func NewTaxonomyFromNCBI(file string) (*Taxonomy, error) {
	return NewTaxonomy(file, 1, 3)
}

// NewTaxonomy only loads nodes from nodes.dmp file.
func NewTaxonomy(file string, childColumn int, parentColumn int) (*Taxonomy, error) {
	if childColumn < 1 || parentColumn < 1 {
		return nil, ErrIllegalColumnIndex
	}

	maxColumns := maxInt(childColumn, parentColumn)

	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, fmt.Errorf("taxdump: %s", err)
	}
	defer func() {
		fh.Close()
	}()

	nodes := make(map[uint32]uint32, 1024)

	n := maxColumns + 1

	childColumn--
	parentColumn--

	items := make([]string, n)
	scanner := bufio.NewScanner(fh)
	var _child, _parent int
	var child, parent uint32
	var maxTaxid uint32
	var root uint32
	for scanner.Scan() {
		stringSplitN(scanner.Text(), "\t", n, &items)
		if len(items) < maxColumns {
			continue
		}

		_child, err = strconv.Atoi(items[childColumn])
		if err != nil {
			continue
		}

		_parent, err = strconv.Atoi(items[parentColumn])
		if err != nil {
			continue
		}

		child, parent = uint32(_child), uint32(_parent)

		// ----------------------------------

		nodes[child] = parent

		if child == parent {
			root = child
		}
		if child > maxTaxid {
			maxTaxid = child
		}

	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("taxdump: %s", err)
	}

	return &Taxonomy{file: file, Nodes: nodes, rootNode: root, maxTaxid: maxTaxid}, nil
}

// NewTaxonomyWithRankFromNCBI parses Taxonomy from nodes.dmp
// from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz .
func NewTaxonomyWithRankFromNCBI(file string) (*Taxonomy, error) {
	return NewTaxonomyWithRank(file, 1, 3, 5)
}

// NewTaxonomyWithRank loads nodes and ranks from nodes.dmp file.
func NewTaxonomyWithRank(file string, childColumn int, parentColumn int, rankColumn int) (*Taxonomy, error) {
	if childColumn < 1 || parentColumn < 1 || rankColumn < 1 {
		return nil, ErrIllegalColumnIndex
	}

	maxColumns := maxInt(childColumn, parentColumn, rankColumn)

	taxid2rankid := make(map[uint32]uint8, 1024)
	ranks := make([]string, 0, 128)
	rank2rankid := make(map[string]int, 128)
	ranksMap := make(map[string]interface{}, 128)

	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, fmt.Errorf("taxdump: %s", err)
	}
	defer func() {
		fh.Close()
	}()

	nodes := make(map[uint32]uint32, 1024)

	n := maxColumns + 1

	childColumn--
	parentColumn--
	rankColumn--

	items := make([]string, n)
	scanner := bufio.NewScanner(fh)
	var _child, _parent int
	var child, parent uint32
	var maxTaxid uint32
	var rank string
	var ok bool
	var rankid int
	var root uint32
	for scanner.Scan() {
		stringSplitN(scanner.Text(), "\t", n, &items)
		if len(items) < maxColumns {
			continue
		}

		_child, err = strconv.Atoi(items[childColumn])
		if err != nil {
			continue
		}

		_parent, err = strconv.Atoi(items[parentColumn])
		if err != nil {
			continue
		}

		child, parent, rank = uint32(_child), uint32(_parent), items[rankColumn]

		// ----------------------------------

		nodes[child] = parent

		if child == parent {
			root = child
		}
		if child > maxTaxid {
			maxTaxid = child
		}

		if rankid, ok = rank2rankid[rank]; ok {
			taxid2rankid[child] = uint8(rankid)
		} else {
			ranks = append(ranks, rank)
			if len(ranks) > 255 {
				return nil, ErrTooManyRanks
			}
			rank2rankid[rank] = len(ranks) - 1
			taxid2rankid[child] = uint8(len(ranks) - 1)
			ranksMap[rank] = struct{}{}
		}

	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("taxdump: %s", err)
	}

	return &Taxonomy{file: file, Nodes: nodes, rootNode: root, maxTaxid: maxTaxid,
		taxid2rankid: taxid2rankid, ranks: ranks, hasRanks: true, Ranks: ranksMap}, nil
}

// TaxId checks if a TaxId is valid in the database.
// If being merged, the new TaxId will be returned.
// If the TaxId is not found or deleted, false will be returned.
func (t *Taxonomy) TaxId(taxid uint32) (uint32, bool) {
	_, ok := t.Nodes[taxid]
	if ok {
		return taxid, true
	}

	// check if it was merged
	var newtaxid uint32
	if newtaxid, ok = t.MergeNodes[taxid]; ok {
		return newtaxid, true
	}

	// // taxid not found, check if it was deleted
	// if _, ok = t.DelNodes[taxid]; ok {
	// 	return taxid, false
	// }

	// not found
	return taxid, false
}

// Name returns the name of a TaxId.
// If being merged, the name of the new TaxId will be returned.
// If the TaxId is not found or deleted, empty will be returned.
func (t *Taxonomy) Name(taxid uint32) string {
	name, ok := t.Names[taxid]
	if ok {
		return name
	}

	// check if it was merged
	var newtaxid uint32
	if newtaxid, ok = t.MergeNodes[taxid]; ok {
		return t.Names[newtaxid]
	}

	return ""
}

// Rank returns the rank of a taxid.
// If being merged, the rank of the new TaxId will be returned.
// If the TaxId is not found or deleted, empty will be returned.
func (t *Taxonomy) Rank(taxid uint32) string {
	if !t.hasRanks {
		panic(ErrRankNotLoaded)
	}
	if i, ok := t.taxid2rankid[taxid]; ok {
		return t.ranks[int(i)]
	}

	if newTaxid, ok := t.MergeNodes[taxid]; ok { // merged
		if i, ok := t.taxid2rankid[newTaxid]; ok {
			return t.ranks[int(i)]
		}
	}

	return "" // taxid not found or deleted
}

// AtOrBelowRank returns whether a taxid is at or below one rank.
func (t *Taxonomy) AtOrBelowRank(taxid uint32, rank string) bool {
	if !t.hasRanks {
		panic(ErrRankNotLoaded)
	}
	var ok bool
	var i uint8

	rank = strings.ToLower(rank)
	if _, ok = t.Ranks[rank]; !ok {
		return false
	}

	if i, ok = t.taxid2rankid[taxid]; ok {
		if rank == t.ranks[int(i)] {
			return true
		}
	}

	// continue searching towards to root node
	var child, parent, newtaxid uint32

	child = taxid
	for {
		parent, ok = t.Nodes[child]
		if !ok { // taxid not found
			// // check if it was deleted
			// if _, ok = t.DelNodes[child]; ok {
			// 	return false
			// }

			// check if it was merged
			if newtaxid, ok = t.MergeNodes[child]; ok {
				child = newtaxid

				if rank == t.ranks[t.taxid2rankid[child]] {
					return true
				}

				parent = t.Nodes[child]
			} else { // not found
				return false
			}
		}

		if parent == 1 {
			break
		}

		if rank == t.ranks[t.taxid2rankid[parent]] {
			return true
		}

		child = parent
	}

	return false
}

// LoadNamesFromNCBI loads scientific names from NCBI names.dmp
func (t *Taxonomy) LoadNamesFromNCBI(file string) error {
	return t.LoadNames(file, 1, 3, 7, "scientific name")
}

// LoadNames loads names.
func (t *Taxonomy) LoadNames(file string, taxidColumn int, nameColumn int, typeColumn int, _type string) error {
	if taxidColumn < 1 || nameColumn < 1 || typeColumn < 1 {
		return ErrIllegalColumnIndex
	}

	maxColumns := maxInt(nameColumn, nameColumn, typeColumn)

	fh, err := xopen.Ropen(file)
	if err != nil {
		return fmt.Errorf("taxdump: %s", err)
	}
	defer func() {
		fh.Close()
	}()

	m := make(map[uint32]string, 1024)

	n := maxColumns + 1

	taxidColumn--
	nameColumn--
	typeColumn--

	filterByType := _type != ""

	items := make([]string, n)
	scanner := bufio.NewScanner(fh)
	var taxid uint64
	for scanner.Scan() {
		stringSplitN(scanner.Text(), "\t", n, &items)
		if len(items) < maxColumns {
			continue
		}

		if filterByType && items[typeColumn] != _type {
			continue
		}

		taxid, err = strconv.ParseUint(items[taxidColumn], 10, 32)
		if err != nil {
			continue
		}

		m[uint32(taxid)] = items[nameColumn]
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("taxdump: %s", err)
	}

	t.Names = m
	t.hasNames = true
	return nil
}

// LoadMergedNodesFromNCBI loads merged nodes from  NCBI merged.dmp.
func (t *Taxonomy) LoadMergedNodesFromNCBI(file string) error {
	return t.LoadMergedNodes(file, 1, 3)
}

// LoadMergedNodes loads merged nodes.
func (t *Taxonomy) LoadMergedNodes(file string, oldColumn int, newColumn int) error {
	if oldColumn < 1 || newColumn < 1 {
		return ErrIllegalColumnIndex
	}

	maxColumns := maxInt(oldColumn, newColumn)

	fh, err := xopen.Ropen(file)
	if err != nil {
		if err == xopen.ErrNoContent {
			return nil
		}
		return fmt.Errorf("taxdump: %s", err)
	}
	defer func() {
		fh.Close()
	}()

	m := make(map[uint32]uint32, 1024)

	n := maxColumns + 1

	oldColumn--
	newColumn--

	items := make([]string, n)
	scanner := bufio.NewScanner(fh)
	var from, to int
	for scanner.Scan() {
		stringSplitN(scanner.Text(), "\t", n, &items)
		if len(items) < maxColumns {
			continue
		}
		from, err = strconv.Atoi(items[oldColumn])
		if err != nil {
			continue
		}
		to, err = strconv.Atoi(items[newColumn])
		if err != nil {
			continue
		}

		m[uint32(from)] = uint32(to)
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("taxdump: %s", err)
	}

	t.MergeNodes = m
	t.hasMergeNodes = true
	return nil
}

// LoadDeletedNodesFromNCBI loads deleted nodes from NCBI delnodes.dmp.
func (t *Taxonomy) LoadDeletedNodesFromNCBI(file string) error {
	return t.LoadDeletedNodes(file, 1)
}

// LoadDeletedNodes loads deleted nodes.
func (t *Taxonomy) LoadDeletedNodes(file string, column int) error {
	if column < 1 {
		return ErrIllegalColumnIndex
	}

	fh, err := xopen.Ropen(file)
	if err != nil {
		if err == xopen.ErrNoContent {
			return nil
		}
		return fmt.Errorf("taxdump: %s", err)
	}
	defer func() {
		fh.Close()
	}()

	m := make(map[uint32]struct{}, 1024)

	maxColumns := column
	n := maxColumns + 1

	column--
	items := make([]string, n)
	scanner := bufio.NewScanner(fh)
	var id int
	for scanner.Scan() {
		stringSplitN(scanner.Text(), "\t", n, &items)
		if len(items) < maxColumns {
			continue
		}
		id, err = strconv.Atoi(items[column])
		if err != nil {
			continue
		}

		m[uint32(id)] = struct{}{}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("taxdump: %s", err)
	}

	t.DelNodes = m
	t.hasDelNodes = true
	return nil
}

// MaxTaxid returns maximum taxid
func (t *Taxonomy) MaxTaxid() uint32 {
	return t.maxTaxid
}

// CacheLCA tells to cache every LCA query result
func (t *Taxonomy) CacheLCA() {
	t.cacheLCA = true
}

// LCA returns the Lowest Common Ancestor of two nodes, 0 for unknown taxid.
func (t *Taxonomy) LCA(a uint32, b uint32) uint32 {
	if a == 0 || b == 0 {
		return 0
	}
	if a == b {
		return a
	}

	// check cache
	var ok bool

	var query uint64
	var tmp interface{}
	if t.cacheLCA {
		query = pack2uint32(a, b)

		tmp, ok = t.lcaCache.Load(query)
		if ok {
			return tmp.(uint32)
		}
	}

	mA := make(map[uint32]struct{}, 16)

	var child, parent, newTaxid uint32
	var flag bool

	child = a
	for {
		parent, ok = t.Nodes[child]
		if !ok {
			flag = false
			if newTaxid, ok = t.MergeNodes[child]; ok { // merged
				child = newTaxid // update child

				parent, ok = t.Nodes[child]
				if ok {
					flag = true
				}
			}

			if !flag {
				if t.cacheLCA {
					t.lcaCache.Store(query, uint32(0))
				}
				return 0
			}
		}
		if parent == child { // root
			mA[parent] = struct{}{}
			break
		}
		if parent == b { // b is ancestor of a
			if t.cacheLCA {
				t.lcaCache.Store(query, b)
			}
			return b
		}
		mA[parent] = struct{}{}

		child = parent
	}

	child = b
	for {
		parent, ok = t.Nodes[child]
		if !ok {
			flag = false
			if newTaxid, ok = t.MergeNodes[child]; ok { // merged
				child = newTaxid // update child

				parent, ok = t.Nodes[child]
				if ok {
					flag = true
				}
			}

			if !flag {
				if t.cacheLCA {
					t.lcaCache.Store(query, uint32(0))
				}
				return 0
			}
		}

		if parent == child { // root
			break
		}
		if parent == a { // a is ancestor of b
			if t.cacheLCA {
				t.lcaCache.Store(query, a)
			}
			return a
		}
		if _, ok = mA[parent]; ok {
			if t.cacheLCA {
				t.lcaCache.Store(query, parent)
			}
			return parent
		}

		child = parent
	}
	return t.rootNode
}

// LineageNames returns nodes' names of the the complete lineage.
func (t *Taxonomy) LineageNames(taxid uint32) []string {
	taxids := t.LineageTaxIds(taxid)
	if taxids == nil {
		return nil
	}

	if !t.hasNames {
		panic(ErrNamesNotLoaded)
	}

	names := make([]string, len(taxids))
	for i, tax := range taxids {
		names[i] = t.Names[tax]
	}
	return names
}

// LineageTaxIds returns nodes' taxid of the the complete lineage.
func (t *Taxonomy) LineageTaxIds(taxid uint32) []uint32 {
	var child, parent, newtaxid uint32
	var ok bool

	child = taxid
	list := make([]uint32, 0, 16)
	for {
		parent, ok = t.Nodes[child]
		if !ok { // taxid not found
			// // check if it was deleted
			// if _, ok = t.DelNodes[child]; ok {
			// 	return nil
			// }

			// check if it was merged
			if newtaxid, ok = t.MergeNodes[child]; ok {
				child = newtaxid
				parent = t.Nodes[child]
			} else { // not found
				return nil
			}
		}

		list = append(list, child)

		if parent == 1 {
			break
		}
		child = parent
	}

	// reversing
	for i, j := 0, len(list)-1; i < j; i, j = i+1, j-1 {
		list[i], list[j] = list[j], list[i]
	}

	return list
}

func pack2uint32(a uint32, b uint32) uint64 {
	if a < b {
		return (uint64(a) << 32) | uint64(b)
	}
	return (uint64(b) << 32) | uint64(a)
}

func minInt(a int, vals ...int) int {
	min := a
	for _, v := range vals {
		if v < min {
			min = v
		}
	}
	return min
}

func maxInt(a int, vals ...int) int {
	min := a
	for _, v := range vals {
		if v > min {
			min = v
		}
	}
	return min
}

func stringSplitN(s string, sep string, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.Index(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+len(sep):]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}
