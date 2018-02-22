#!/usr/bin/env bash
set -e

docker exec -ti ambrosejcarr/seqc:latest nose2 seqc.test.TestSEQC.test_local

