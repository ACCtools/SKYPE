import logging

# 파이프라인 전역 logging 수준. 모든 SKYPE 스크립트가 'from skype_utils import *' 를 가장 먼저
# 실행하므로, 여기서 basicConfig 를 호출해두면 각 스크립트의 자체 basicConfig 보다 먼저 root
# logger 가 설정되어 이 레벨이 전체에 적용된다(자체 basicConfig 는 핸들러가 이미 있어 무시됨).
# 디버그 출력(logging.debug)을 보려면 logging.DEBUG 로 바꾼다.
LOG_LEVEL = logging.INFO
logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=LOG_LEVEL,
    datefmt='%m/%d/%Y %I:%M:%S %p',
    force=True
)

K = 1000
M = 1000 * K

HARD_PATH_COUNT_BASELINE = 500 * K
chrY_MINIMUM_RATIO = 0.7

CHROMOSOME_COUNT = 23
