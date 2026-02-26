#!/usr/bin/env python3

__version__ = '0.4.0'

from argparse import ArgumentParser
from dataclasses import dataclass, field
import glob
import json
from multiprocessing import Pool
import os
import pickle
import subprocess
import sys
import tempfile


DEFAULT_THREADS = 1
DEFAULT_MIN_DUPLEX_DEPTH = 2
DEFAULT_MIN_BASE_QUALITY = 30

DIR_COV = 'cov'
DIR_DSA = 'dsa'
DIR_PART = 'part'

DONE_FILE_NAME = 'done'
BED_FILE_NAME = 'ranges.bed'
DSA_FILE_NAME = 'dsa.bed.gz'
DSA_REPORT_FILE_NAME = 'report.json'


class GInterval:
    def __init__(self, chrr: str, beg: int, end: int):
        self.chr = chrr
        if (end < beg):
            raise ValueError("Interval %s: %s - %s is invalid!" %
                             (chrr, beg, end))
        self.beg = beg
        self.end = end
        self.l = end - beg + 1

    def convert2DSAInput(self):
        # zero based and inclusive of end
        return(GInterval(self.chr, self.beg - 1, self.end - 1))

    def write_bed(self, fh):
        r = self.convert2DSAInput()
        fh.write(f"{r.chr}\t{r.beg}\t{r.end}\n")


def get_part_file(tmp_dir: str, fn: str) -> str:
    return os.path.join(tmp_dir, DIR_PART, fn)


def get_cov_file(tmp_dir: str, fn: str) -> str:
    return os.path.join(tmp_dir, DIR_COV, fn)


def get_dsa_file(tmp_dir: str, fn: str) -> str:
    return os.path.join(tmp_dir, DIR_DSA, fn)


@dataclass(slots=True)
class Cmd:
    array_index: int
    cmd: str

    def __post_init__(self) -> None:
        assert self.array_index >= 0
        assert self.cmd

    def run(self, dry: bool = False) -> None:
        if not dry:
            p = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            _, stderr = p.communicate()
            if p.returncode != 0:
                error = stderr.decode()
                sys.stderr.write("Error processing at array index %d!\n" % self.array_index)
                raise ValueError(error)
        else:
            print(self.cmd)

    def push_and(self, cmd: str):
        self.cmd += ' && ' + cmd


def chain_cmds(cmds: list[str]) -> str:
    return ' && '.join(cmds)


@dataclass(slots=True, frozen=True)
class JobInfo:
    dir: str
    array_index: int

    @property
    def work_dir(self) -> str:
        # e.g., /tmp/dsa/3/
        return os.path.join(self.dir, str(self.array_index))

    def get_fp(self, fn: str) -> str:
        return os.path.join(self.dir, fn)

    def get_work_dir_fp(self, fn: str) -> str:
        return os.path.join(self.work_dir, fn)

    def _qualify_fn(self, fn: str) -> str:
        return f"{self.array_index}.{fn}"

    @property
    def dsa_fp(self) -> str:
        # e.g., 1/dsa.bed.gz
        return self.get_work_dir_fp(DSA_FILE_NAME)

    @property
    def dsa_ln(self) -> str:
        # e.g., 1.dsa.bed.gz -> 1/dsa.bed.gz
        return self.get_work_dir_fp(self._qualify_fn(DSA_FILE_NAME))

    @property
    def done_fp(self) -> str:
        # e.g., 1.done
        return self.get_work_dir_fp(self._qualify_fn(DONE_FILE_NAME))

    @property
    def ranges_fp(self) -> str:
        return os.path.join(self.work_dir, BED_FILE_NAME)


@dataclass(slots=True)
class CmdBuilder:
    exe: str
    job: JobInfo

    _tokens: list[str] = field(init=False)

    def push(self, a: str) -> None:
        self._tokens.append(a)

    def push_option(self, k: str, v: str | int) -> None:
        self.push(f"-{k}")
        # TODO: consider quoting strings?
        self.push(str(v))

    def cmd(self) -> Cmd:
        s = ' '.join([self.exe, *self._tokens])
        cmd = Cmd(self.job.array_index, s)

        return cmd


@dataclass(slots=True, frozen=True)
class DSAArgs:
    exe: str
    normal: str
    duplex: str
    ref: str
    q: int
    d: int
    mapQ: int
    no_test: bool
    snp_mask: str | None
    noise_mask: str | None


def main():
    # TODO: generate BED file (in the output directory for ease of QC?)
    pass


def get_symlink_cmd(src: str, dest: str) -> str:
    return f"ln -s {src} {dest}"


def prepare_dsa_job(args: DSAArgs, job: JobInfo) -> Cmd | None:
    "Generates the command and prepares the working directory"

    # check for restarts
    if (
        os.path.isfile(job.dsa_ln) and
        os.path.isfile(job.done_fp)
    ):
        return None

    ranges_fp = job.get_work_dir_fp(BED_FILE_NAME)
    with open(ranges_fp, 'w') as fh:
        for interval in intervalsPerCPU[i]:
            interval.write_bed(fh)

    b = CmdBuilder(args.exe, job)
    b.push_option('A', args.normal)
    b.push_option('B', args.duplex)
    b.push_option('I', job.ranges_fp)
    b.push_option('R', args.ref)
    b.push_option('d', args.d)
    b.push_option('Q', args.q)
    b.push_option('M', args.mapQ)

    if args.no_test:
        b.push('-t')

    if args.snp_mask:
        b.push_option('C', args.snp_mask)

    if args.noise_mask:
        b.push_option('D', args.noise_mask)

    cmd = b.cmd()

    # Create symlink for the dsa.bed.gz file as <INDEX>.dsa.bed.gz in the dsa directory
    cmd.push_and(get_symlink_cmd(job.dsa_fp, job.dsa_ln))

    # Create the done file in the dsa directory
    cmd.push_and(f"touch {job.done_fp}")
    return cmd


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-j', '--index', type=int, help='index of the LSF job array. One based')
    p.add_argument('-k', '--max_index', type=int, help='maximum index of the LSF job array')
    p.add_argument('-t', '--threads', type=int, default=DEFAULT_THREADS, help=f"number of threads ({DEFAULT_THREADS})")

    p.add_argument('-R', '--ref', required=True, help="referene sequence")
    p.add_argument('-A', '--normal', required=True, help="normal BAM / CRAM")
    p.add_argument('-B', '--duplex', required=True, help="duplex (tumour) BAM / CRAM")

    p.add_argument('-C', '--snp', help="SNP BED (gz) file")
    p.add_argument('-D', '--mask', help="mask BED (gz) file")
    p.add_argument('-d', type=int, default=DEFAULT_MIN_DUPLEX_DEPTH, help=f"minimum duplex depth ({DEFAULT_MIN_DUPLEX_DEPTH})")
    p.add_argument('-q', type=int, default=DEFAULT_MIN_BASE_QUALITY, help=f"minimum base quality for normal ({DEFAULT_MIN_BASE_QUALITY}])")
    p.add_argument('--no_test', action='store_true', help="skip BAM format tests, use with caution")
    p.add_argument('--dry', action='store_true', help="print the commands and exit")

    p.add_argument('--out', default='.', help='path of the output files and scratch directory (.)')
    p.add_argument('-v', '--version', action='version', version=__version__)
    args = p.parse_args()

    exe: str | None = os.getenv('DSA_EXE')
    assert exe

    if not os.path.isdir(args.out):
        p.error("Specified out directory %s is not accessible!" % args.out)

    try:
        testfile = tempfile.TemporaryFile(dir=args.out)
        testfile.close()
    except OSError:
        sys.exit("\nCan't write to out directory %s\n" % args.out)

    tmp_dir = os.path.join(args.out, 'tmpNanoSeq')

    # Validate preceding 'part' step
    part_args_fp = get_part_file(tmp_dir, 'args.json')
    if not os.path.isfile(part_args_fp):
        sys.exit("\nMust run cov and part submmands prior to dsa\n")
    else:
        with open(part_args_fp) as iofile:
            oargs = json.load(iofile)
    njobs: int = oargs['jobs']

    assert isinstance(njobs, int)

    part_job = JobInfo(os.path.join(tmp_dir, DIR_PART), 1)
    part_done_fp = part_job.done_fp
    # part_done_fp = get_part_file(tmp_dir, '1.done')
    # part_ipc_fp = get_part_file(tmp_dir, 'intervalsPerCPU.dat')
    part_ipc_fp = part_job.get_fp('intervalsPerCPU.dat')

    if len(glob.glob(part_done_fp)) != 1:
        sys.exit("\npart job did not complete correctly\n")
    if len(glob.glob(part_ipc_fp)) != 1:
        sys.exit("\npart job did not complete correctly\n")

    # make sure that number of jobs matches what was specified in part
    if (args.max_index is not None):
        # array execution
        if (args.max_index < njobs):
            sys.exit(
                "\nLSF array size must match number of jobs specified for part (%s)\n" % njobs)
        if (args.index > njobs):
            print(
                "\nWarning specified LSF array size is larger than jobs specified for part (%s)\n" % njobs)
            sys.exit(0)
    else:
        # multithread
        if (args.threads < njobs):
            sys.exit(
                "\nNumber of threads must match number of jobs specified for part (%s)\n" % njobs)
        if (args.threads > njobs):
            print(
                "\nWarning number of threads is larger than jobs specified for part (%s)\n" % njobs)

    with open(part_ipc_fp, 'rb') as iofile:
        intervalsPerCPU = pickle.load(iofile)

    mapQ = None
    cov_job = JobInfo(os.path.join(tmp_dir, DIR_COV), -1)
    cov_args_fp = cov_job.get_fp('args.json')
    with open(cov_args_fp) as iofile:
        mapQ = json.load(iofile)['Q']

    # TODO: generate!
    ranges_fp = ""

    a = DSAArgs(
        exe=exe,
        normal=args.normal,
        duplex=args.duplex,
        ref=args.ref,
        d=args.d,
        q=args.q,
        mapQ=mapQ,
        snp_mask=args.snp or None,
        noise_mask=args.mask or None,
        no_test=args.no_test)

    # TODO: moving the dsa table to the root could facilitate clean-up...? Unless it's done at the directory level.

    if args.index is None or args.index == 1:
        with open("%s/dsa/nfiles" % (tmp_dir), "w") as iofile:
            iofile.write(str(njobs))

    # execute dsa commans
    print("Starting dsa calculation\n")
    dsa_dir = os.path.join(tmp_dir, DIR_DSA)
    if (args.index is None):
        commands: list[Cmd] = []
        for i in range(njobs):
            cmd = prepare_dsa_job(a, JobInfo(dsa_dir, i + 1))
            if cmd is not None:
                commands.append(cmd)

        # multithread
        with Pool(args.threads) as p:
            p.map(lambda x: x.run(dry=args.dry), commands)

    else:
        # array execution
        cmd = prepare_dsa_job(a, JobInfo(dsa_dir, args.index + 1))
        if cmd:
            cmd.run(dry=args.dry)

    print("Completed dsa calculation\n")
