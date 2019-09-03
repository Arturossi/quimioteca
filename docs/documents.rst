.. Install inscrutions for 'quimioteca'

===========================================
Instruções de instalação do site QUIMIOTECA
===========================================

Neste guia serão listados os passos necessários para se instalar completamente o servidor da quimioteca brasileira.

Pré-requisitos
==============

* Sistema Operacional (apenas um)
    * Ubuntu 18.04
    * CentOS 7

* Bibliotecas necessárias
    * git
    * Conda ou Miniconda
    * postgresql
    * postgresql-contrib

Procedimentos
=============

Aqui serão listados os procedimentos para instalação do sistema, contendo códigos, partindo de uma instalação limpa do sistema Ubuntu 18.04 / cent-os rhel fedora.

No Ubuntu
---------

Primeiramente deve-se atualizar as informações dos pacotes

.. code-block:: console

    sudo apt update

Então se instala o postgreesql e o git

.. code-block:: console

    sudo apt install -y postgresql postgresql-contrib libxrender-dev


No CentOS
---------

Primeiramente deve-se atualizar as informações dos pacotes

.. code-block:: console

    sudo yum check-update

Então se instala o postgreesql e o git
.. code-block:: console

    sudo yum install -y postgresql-server postgresql-contrib git libxrender-dev


Esses passos são comuns para ambos sistemas
-------------------------------------------

Baixe o Miniconda (atenção para o link de download, use o mais recente)

.. code-block:: console

    wget -P /tmp/ https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

Execute o instalador

.. code-block:: console

    sudo chmod 755 /tmp/Miniconda3-latest-Linux-x86_64.sh && /tmp/Miniconda3-latest-Linux-x86_64.sh

Crie a base de dados, o usuário e dê privilégios a ele na base criada (a criação é feita pelo usuário que o postgree cria "postgres")

.. DANGER::
    Quando usar os comandos do postgree sempre se lembre de colocar um ; (Ponto e vírgula) no final, ou então você terá problemas!

Use o comando psql através do usuário postgres

.. code-block:: console

    sudo -u postgres psql

Já logado como o usuário postgres crie a base de dados

.. code-block:: console

    create database quimioteca;

Crie o usuário

.. code-block:: console

    create user quimioteca with encrypted password '123456';

Dê os privilégios administrativos para o usuário criado

.. code-block:: console

    grant all privileges on database quimioteca to quimioteca;

Saia do ambiente do postgree

.. code-block:: console

    \q

Crie uma pasta para o site no home

.. code-block:: console

    mkdir ~/sites/ && cd ~/sites

Agora baixe do git a última versão do site

.. code-block:: console

    git clone https://github.com/Arturossi/quimioteca.git

Instale o ambiente virtual através do arquivo environment.yml

.. code-block:: console

    conda env create -f environment.yml

Ative o ambiente criado

.. code-block:: console

    conda activate chemdb


Para criar acesso ao servidor externamente
------------------------------------------

Vamos descobrir qual IP externo que temos (grave esse IP)

.. code-block:: console

    dig +short myip.opendns.com @resolver1.opendns.com

Entre no diretório do site

.. code-block:: console

    cd quimioteca

Vamos permitir o site no IP acima
.. code-block:: console

    vim ./quimioteca/settings.py

Por volta da linha 28 tem-se o seguinte comando

.. code-block:: console

    ALLOWED_HOSTS = []

Dentro do [] adicione o seu ip na forma de string. Ex: suponha que seu IP seja 123.45.67.890

.. code-block:: console

    ALLOWED_HOSTS = ['123.45.67.890']

Agora vamos migrar o banco de dados (caso algum erro ocorra, renomeie ou remova (eu prefiro renomear para migrations.bak por segurança) a pasta ./chemo/migrations

.. code-block:: console

    ./manage.py makemigrations

Se tudo sair ok, dê o comando

.. code-block:: console

    ./manage.py migrate

Agora você já criou as tabelas no banco, vamos iniciar o servidor na porta 8000

.. code-block:: console

    ./manage.py runserver 0:8000

Caso você esteja acessando o servidor via ssh e queira que o servidor continue rodando mesmo após finalizar a sessão ssh, use este comando abaixo ao invés do acima

.. code-block:: console

    nohup ./manage.py runserver 0:8000 &

