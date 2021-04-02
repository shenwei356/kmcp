// Copyright © 2020-2021 Wei Shen <shenwei356@gmail.com>
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

package cmd

import (
	"fmt"
	"path/filepath"
	"sync"

	"github.com/shenwei356/unikmer"
	"github.com/shenwei356/util/pathutil"
)

func loadTaxonomy(opt *Options, path string) *unikmer.Taxonomy {
	if opt.Verbose {
		log.Infof("loading Taxonomy from: %s", path)
	}
	var t *unikmer.Taxonomy
	var err error

	t, err = unikmer.NewTaxonomyWithRankFromNCBI(filepath.Join(path, "nodes.dmp"))
	if err != nil {
		checkError(fmt.Errorf("err on loading Taxonomy nodes: %s", err))
	}

	if opt.Verbose {
		log.Infof("  %d nodes in %d ranks loaded", len(t.Nodes), len(t.Ranks))
	}

	var existed bool

	var wg sync.WaitGroup
	wg.Add(3)

	go func() {
		defer wg.Done()
		existed, err = pathutil.Exists(filepath.Join(path, "names.dmp"))
		if err != nil {
			checkError(fmt.Errorf("err on checking file names.dmp: %s", err))
		}
		if existed {
			err = t.LoadNamesFromNCBI(filepath.Join(path, "names.dmp"))
			if err != nil {
				checkError(fmt.Errorf("err on loading Taxonomy names: %s", err))
			}
		} else {
			checkError(fmt.Errorf("names.dmp not found: %s", err))
		}
		if opt.Verbose {
			log.Infof("  %d names loaded", len(t.Names))
		}
	}()

	go func() {
		defer wg.Done()
		existed, err = pathutil.Exists(filepath.Join(path, "delnodes.dmp"))
		if err != nil {
			checkError(fmt.Errorf("err on checking file merged.dmp: %s", err))
		}
		if existed {
			err = t.LoadDeletedNodesFromNCBI(filepath.Join(path, "delnodes.dmp"))
			if err != nil {
				checkError(fmt.Errorf("err on loading Taxonomy nodes: %s", err))
			}
		}
		if opt.Verbose {
			log.Infof("  %d deleted nodes loaded", len(t.DelNodes))
		}
	}()

	go func() {
		defer wg.Done()
		existed, err = pathutil.Exists(filepath.Join(path, "merged.dmp"))
		if err != nil {
			checkError(fmt.Errorf("err on checking file merged.dmp: %s", err))
		}
		if existed {
			err = t.LoadMergedNodesFromNCBI(filepath.Join(path, "merged.dmp"))
			if err != nil {
				checkError(fmt.Errorf("err on loading Taxonomy merged nodes: %s", err))
			}
		}
		if opt.Verbose {
			log.Infof("  %d merged nodes loaded", len(t.MergeNodes))
		}
	}()

	wg.Wait()

	t.CacheLCA()

	return t
}
