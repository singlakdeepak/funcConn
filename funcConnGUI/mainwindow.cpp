#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "editlist.h"
#include "editfeat.h"
#include "choosecombinations.h"
#include <QSpinBox>
#include <QMessageBox>
#include<QFileDialog>



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->radioButton_Unprocessed->setChecked(true);
    ui->checkBox->setChecked(false);
    ui->checkBox_2->setChecked(false);
    ui->checkBox_3->setChecked(false);
    ui->spinBox->setRange(1,5);
    ui->radioButton_wt_Groups->setChecked(true);
    ui->chooseCombsButton->hide();

    QString DefAnalysisName = "FConnectivityAnl";
    QString FSLDIRpos1 = "/usr/share/fsl/5.0/bin/";
    QString FSLDIRpos2 = "/usr/share/fsl/bin/";

    if (QDir(FSLDIRpos1).exists()){
        ui->lineEdit_FSLDIR->setText(FSLDIRpos1);
    }
    else if (QDir(FSLDIRpos2).exists()){
        ui->lineEdit_FSLDIR->setText(FSLDIRpos2);
    }

    ui->lineEdit_AnalysisName->setText(DefAnalysisName);
    ui->lineEdit_OutDir->setText(OutDir + DefAnalysisName);
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_chooseData_clicked()
{
    QString spin = ui->spinBox->text();
    int spin_1 = spin.toInt();
    if (ui->radioButton_PreprocwtFSL->isChecked()){
        editFeat editfeat;
        editfeat.setGroups(spin_1);
        editfeat.setModal(true);
        editfeat.exec();
    }
    else {
        EditList editlist;
        editlist.setGroups(spin_1);
        editlist.setModal(true);
        editlist.exec();
    }

}

void MainWindow::on_pushButton_7_clicked()
{
    ui->statusBar->showMessage("Put the reference MNI structural file here.",10000); //It will show message for 10 s.
}

void MainWindow::on_pushButton_3_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_6_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_12_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_18_clicked()
{
    QMessageBox::information(this,"Help", "If correlation files given, then there is no requirement of input.");
}

void MainWindow::on_commandLinkButton_Reference_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_Reference->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
            new_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
    ui->lineEdit_Reference->setText(new_input);

    }
    else {
        ui->lineEdit_Reference->setText(already_input);
    }
}

void MainWindow::on_commandLinkButton_OutDir_clicked()
{
    QString new_Dir;

    if (QDir(OutDir).exists())
    {
        new_Dir=QFileDialog::getExistingDirectory(
                                                    this,
                                                    tr("Choose Directory"),
                                                    OutDir);
    }
    else {
        new_Dir =QFileDialog::getExistingDirectory(
                                                this,
                                                tr("Choose Directory"),
                                                CommonDir);
    }
    if (new_Dir!=NULL){
        ui->lineEdit_OutDir->setText(new_Dir +"/"+ ui->lineEdit_AnalysisName->text());
        OutDir = new_Dir+"/";
    }
    else{
        ui->lineEdit_OutDir->setText(OutDir + ui->lineEdit_AnalysisName->text());
    }
}

void MainWindow::on_commandLinkButton_ROIFile_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_ROIFile->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
        already_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        already_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        ui->lineEdit_ROIFile->setText(new_input);

    }
    else {
        ui->lineEdit_ROIFile->setText(already_input);
    }
}

void MainWindow::on_commandLinkButton_FSLDIR_clicked()
{

    QString filename=QFileDialog::getExistingDirectory(
                                                this,
                                                tr("Choose Directory"),
                                                "/home/deepak/");
    ui->lineEdit_FSLDIR->setText(filename);
}



void MainWindow::on_commandLinkButton_CorrFile_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_CorrFile->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
      already_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        already_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        ui->lineEdit_CorrFile->setText(new_input);

    }
    else {
        ui->lineEdit_CorrFile->setText(already_input);
    }
}

void MainWindow::on_help_ROIFile_clicked()
{
    QMessageBox::information(this,"Caution", "This file must be registered to reference.");
}


void MainWindow::on_lineEdit_AnalysisName_textEdited(const QString &arg1)
{
    ui->lineEdit_OutDir->setText(OutDir + arg1);

}

void MainWindow::on_radioButton_bw_Groups_clicked(bool checked)
{

    if (checked){
         ui->chooseCombsButton->show();
     }
}



void MainWindow::on_radioButton_wt_Groups_clicked(bool checked)
{
    if (checked){
        ui->chooseCombsButton->hide();
    }
}

void MainWindow::on_chooseCombsButton_clicked()
{
    chooseCombinations choosecombinations;
    choosecombinations.setModal(true);
    choosecombinations.exec();
}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    if (arg1==1){
        ui->radioButton_bw_Groups->hide();
    }
    else {
        ui->radioButton_bw_Groups->show();
    }
}
