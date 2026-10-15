# GitHub 업로드와 실행

## 업로드 범위

`crtmedsim` 폴더의 내용이 저장소의 최상위에 오도록 업로드한다.
공개 파일은 패키지 소스, 고정된 시뮬레이션 엔진, 테스트, 실행 예제,
설명 문서로 구성된다. 라이선스는 추후 선택할 수 있도록 보류 상태이며,
`DESCRIPTION`의 연락처는 프로젝트용 예시 주소다.

## 로컬 설치와 첫 실행

터미널에서 `crtmedsim` 폴더로 이동한 뒤 실행한다.

```sh
Rscript scripts/install_dependencies.R
R CMD INSTALL .
Rscript scripts/reproduce.R smoke runs/smoke 1
```

첫 실행은 K=50, 반복 1회이고, 방법과 시나리오는 모두 포함한다.
적합 및 적분 설정은 최종 연구 설정을 유지한다.
전체 연구 실행은 다음 명령을 사용한다.

```sh
Rscript scripts/reproduce.R full runs/full 10
```

전체 설정은 K=50, 75, 100, 150 각각 반복 2,000회다.
동일 명령을 다시 실행하면 완료된 결과를 확인하고 남은 반복을 실행한다.
로컬에서 다른 계산과 함께 실행할 때는 마지막 숫자를 `1`로 둔다.

## GitHub 저장소 생성 후 업로드

GitHub에서 빈 저장소를 생성하고, 해당 저장소의 URL을 사용한다.

```sh
git init
git add .
git diff --cached --stat
git commit -m "Add simulation reproduction package"
git branch -M main
git remote add origin <repository-url>
git push -u origin main
```

업로드 전 `git diff --cached --stat`에서 파일 목록을 확인한다.
`.gitignore`가 실행 결과, RDS 파일, 로그, 압축 파일을 제외한다.
Git 커밋의 작성자 정보와 이메일은 Git 설정에서 별도로 관리한다.
GitHub 계정의 이메일 공개 설정과 비공개 커밋 이메일도 확인한다.

## 서버 실행과 결과 전달

서버 배열 실행은 `docs/SLURM.md`의 명령을 사용한다.
완료 후 결과 분석에는 다음 파일 하나를 전달하면 된다.

```text
runs/full/results/annals_simulation_bundle.rds
```

수치 표, 그림, LaTeX 표는 `runs/full/results/report/`에 생성된다.

