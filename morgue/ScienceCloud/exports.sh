export OS_AUTH_URL=https://cloud.s3it.uzh.ch:5000/v2.0
export OS_USERNAME=mbarbo
export OS_TENANT_NAME=bascompte.ieu.mnf.uzh
export OS_PROJECT_NAME=bascompte.ieu.mnf.uzh
read -p "Password: " -s mypassword
export OS_PASSWORD=$mypassword
unset OS_REGION_NAME
export PYTHONPATH=$PWD
