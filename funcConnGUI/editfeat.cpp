#include "editfeat.h"
#include "ui_editfeat.h"
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include<QFileDialog>



editFeat::editFeat(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::editFeat)
{
    ui->setupUi(this);
}

editFeat::~editFeat()
{
    delete ui;
}

void editFeat::setGroups(int groups)
{
    for (int i =1; i<=groups;i++){
        ui->comboBox_Groups->addItem("Group " + QString::number(i));
    }
}

void editFeat::on_pushButton_Save_clicked()
{
    int Groupindex = ui->comboBox_Groups->currentIndex();

    QFile file(WorkingDir+"/Group_"+QString::number(Groupindex +1)+"_FeatDirs.txt");
    if (!file.open(QFile::WriteOnly | QFile::Text)){
        QMessageBox::warning(this, "title","file not open");
    }
    QTextStream out(&file);
    QString text = ui->TextEdit_FeatDirs->toPlainText();
    out << text;
    file.flush();
    file.close();
}

void editFeat::on_comboBox_Groups_currentIndexChanged(int index)
{
    ui->TextEdit_FeatDirs->clear();

    QFile file_Feat( WorkingDir+ "/Group_"+QString::number(index+1)+"_FeatDirs.txt");
    if (file_Feat.open(QFile::ReadOnly | QFile::Text)){
    QTextStream in_Feat(&file_Feat);
    QString text_Feat = in_Feat.readAll();
    ui->TextEdit_FeatDirs->setPlainText(text_Feat);

    file_Feat.close();
    }
}

void editFeat::on_ld_featDirList_clicked()
{
    QFile Fout(currentFeatDir) ;
    QString new_input;
    if (Fout.exists())
    {
        new_input=QFileDialog::getOpenFileName(this,
                                               tr("Open Feat File"),
                                                    currentFeatDir,
                                                    "Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open Feat File"),
                                                CommonDir,
                                                "Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        QFile file(new_input);
        if(!file.open(QIODevice::ReadOnly)) {
            QMessageBox::information(0, "error", file.errorString());
        }
        ui->TextEdit_FeatDirs->clear();
        QTextStream in(&file);

        while(!in.atEnd()) {
            QString line = in.readLine();
            ui->TextEdit_FeatDirs->appendPlainText(line);
        }

        file.close();
        currentFeatDir = new_input;
    }
}
