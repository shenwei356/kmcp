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
	"os"
	"runtime"

	"github.com/spf13/cobra"
)

// RootCmd represents the base command when called without any subcommands
var RootCmd = &cobra.Command{
	Use:   "kmcp",
	Short: "Kmer-based Metagenomics Classification and Profilling",
	Long: fmt.Sprintf(`
    Program: kmcp (K-mer-based Metagenomic Classification and Profiling)
    Version: v%s
  Documents: https://shenwei356.github.io/kmcp
Source code: https://github.com/shenwei356/kmcp

kmcp is a tool for metagenomic classification and profiling.

`, VERSION),
}

// Execute adds all child commands to the root command sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := RootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(-1)
	}
}

func init() {

	defaultThreads := runtime.NumCPU()

	RootCmd.PersistentFlags().IntP("threads", "j", defaultThreads, "number of CPUs to use")
	// RootCmd.PersistentFlags().BoolP("verbose", "", false, "print verbose information (recommended)")
	RootCmd.PersistentFlags().BoolP("quiet", "q", false, "do not print any verbose information")
	RootCmd.PersistentFlags().StringP("infile-list", "i", "", "file of input files list (one file per line), if given, they are appended to files from cli arguments")
}
