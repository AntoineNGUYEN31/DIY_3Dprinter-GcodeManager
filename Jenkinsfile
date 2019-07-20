node('slave-localhost') {
  stage('Prepare'){
    sh "uname -r"
  }      
  stage('Build') {
    input 'Wait'
  }
}
