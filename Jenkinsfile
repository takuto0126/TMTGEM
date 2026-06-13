pipeline {
    agent any

    triggers {
        githubPush()
    }

    stages {

        stage('Checkout') {
            steps {
                checkout scm
            }
        }

        stage('Run tests') {
            steps {
                sh '''
                    chmod +x tests/test_Tohoku_whole.sh
                    ./tests/test_Tohoku_whole.sh
                '''
            }
        }
    }
}