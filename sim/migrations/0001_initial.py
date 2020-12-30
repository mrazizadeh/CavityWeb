# Generated by Django 3.1.4 on 2020-12-30 04:59

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Input',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('A', models.FloatField(default=1.0, verbose_name=' amplitude (m)')),
                ('b', models.FloatField(default=0.0, verbose_name=' damping coefficient (kg/s)')),
                ('w', models.FloatField(default=6.283185307179586, verbose_name=' frequency (1/s)')),
                ('T', models.FloatField(default=18, verbose_name=' time interval (s)')),
            ],
        ),
    ]
