#!/usr/bin/env bash
# Fetch HEK293 comparison summary from EC2 via SSM (or S3 when uploaded).
set -euo pipefail

INSTANCE_ID="${INSTANCE_ID:-i-04818df73a6cbff34}"
WORKDIR="${WORKDIR:-/home/ec2-user/work/hek293_benchmark}"

CMD_ID=$(aws ssm send-command \
  --instance-ids "$INSTANCE_ID" \
  --document-name "AWS-RunShellScript" \
  --parameters "commands=[\"if [ -f ${WORKDIR}/hek293_comparison_summary.txt ]; then cat ${WORKDIR}/hek293_comparison_summary.txt; elif grep -q HEK293_COMPARISON_DONE ${WORKDIR}/run.log 2>/dev/null; then tail -50 ${WORKDIR}/run.log; else echo INCOMPLETE; tail -20 ${WORKDIR}/run.log 2>/dev/null || tail -20 ${WORKDIR}/nohup.out 2>/dev/null; fi\"]" \
  --query 'Command.CommandId' --output text)

sleep 5
aws ssm get-command-invocation \
  --command-id "$CMD_ID" \
  --instance-id "$INSTANCE_ID" \
  --query 'StandardOutputContent' --output text
