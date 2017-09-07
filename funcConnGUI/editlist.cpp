#include "editlist.h"
#include "ui_editlist.h"
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include<QFileDialog>



EditList::EditList(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EditList)
{
    ui->setupUi(this);
    ui->TextEdit_Functional->setPlainText(WorkingDir);
}

EditList::~EditList()
{
    delete ui;
}

void EditList::setGroups(int groups)
{
    for (int i =1; i<=groups;i++){
        ui->comboBox->addItem("Group " + QString::number(i));
    }
}


void EditList::on_pushButton_Save_clicked()
{
    int Groupindex = ui->comboBox->currentIndex();

    if (ui->checkBox->isChecked()){
        bool write = true;
        QString text = ui->TextEdit_Functional->toPlainText();
        QStringList FuncUrls = text.split("\n");
        QString tempFuncString;
        int length_Urls = FuncUrls.size();
        for (int i =0; i < length_Urls; i++){
            tempFuncString = FuncUrls.at(i);
            QFile funcFile(tempFuncString);
            if (!funcFile.exists()){
                QMessageBox::warning(this,"title","File with the following path: ''"+ tempFuncString + "'' doesn't exist. Try again.");
                write = false;
                break;
            }
        }
        if (write){
            QFile file(WorkingDir+"/Group_"+QString::number(Groupindex +1)+"_Functional.txt");
            if (!file.open(QFile::WriteOnly | QFile::Text)){
                QMessageBox::warning(this, "title","file not open");
            }
            QTextStream out(&file);
            out << text;
            file.flush();
            file.close();
        }
    }
    if (ui->checkBox_2->isChecked())
    {
        bool write = true;
        QString text = ui->TextEdit_Structural->toPlainText();
        QStringList StructUrls = text.split("\n");
        QString tempStructString;
        int length_Urls = StructUrls.size();
        for (int i =0; i < length_Urls; i++){
            tempStructString = StructUrls.at(i);
            QFile StructFile(tempStructString);
            if (!StructFile.exists()){
                QMessageBox::warning(this,"title","File with the following path: ''"+ tempStructString + "'' doesn't exist. Try again.");
                write = false;
                break;
            }
        }
        if (write){
            QFile file(WorkingDir+"/Group_"+QString::number(Groupindex +1)+"_Structural.txt");
            if (!file.open(QFile::WriteOnly | QFile::Text)){
                QMessageBox::warning(this, "title","file not open");
            }
            QTextStream out(&file);
            out << text;
            file.flush();
            file.close();
        }
    }
}


void EditList::on_comboBox_currentIndexChanged(int index)
{
    ui->TextEdit_Functional->clear();
    ui->TextEdit_Structural->clear();

    QFile file_Functional( WorkingDir+ "/Group_"+QString::number(index+1)+"_Functional.txt");
    if (file_Functional.open(QFile::ReadOnly | QFile::Text)){


    QTextStream in_Functional(&file_Functional);
    QString text_Functional = in_Functional.readAll();
    ui->TextEdit_Functional->setPlainText(text_Functional);

    file_Functional.close();
    }

    QFile file_Structural(WorkingDir+"/Group_"+QString::number(index+1)+"_Structural.txt");
    if (file_Structural.open(QFile::ReadOnly | QFile::Text)){

    QTextStream in_Structural(&file_Structural);
    QString text_Structural = in_Structural.readAll();
    ui->TextEdit_Structural->setPlainText(text_Structural);

    file_Structural.close();
    }
}

void EditList::on_ld_funcFile_clicked()
{
    QFile Fout(CurrentFuncFile) ;
    QString new_input;
    if (Fout.exists())
    {
        new_input=QFileDialog::getOpenFileName(this,
                                               tr("Open Functional File"),
                                                    CurrentFuncFile,
                                                    "Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open Functional File"),
                                                CommonDir,
                                                "Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        QFile file(new_input);
        if(!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(0, "error", file.errorString());
        }
        ui->TextEdit_Functional->clear();
        QTextStream in(&file);

        while(!in.atEnd()) {
            QString line = in.readLine();
            ui->TextEdit_Functional->appendPlainText(line);
        }

        file.close();
        CurrentFuncFile = new_input;
    }

}

void EditList::on_ld_StrFile_clicked()
{
    QFile Fout(CurrentStructFile) ;
    QString new_input;
    if (Fout.exists())
    {
        new_input=QFileDialog::getOpenFileName(this,
                                               tr("Open Functional File"),
                                                    CurrentStructFile,
                                                    "Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open Functional File"),
                                                CommonDir,
                                                "Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        QFile file(new_input);
        if(!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(0, "error", file.errorString());
        }
        ui->TextEdit_Structural->clear();
        QTextStream in(&file);

        while(!in.atEnd()) {
            QString line = in.readLine();
            ui->TextEdit_Structural->appendPlainText(line);
        }

        file.close();
        CurrentStructFile = new_input;
    }
}
