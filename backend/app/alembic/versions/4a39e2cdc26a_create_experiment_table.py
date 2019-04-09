"""create experiment table

Revision ID: 4a39e2cdc26a
Revises: e6ae69e9dcb9
Create Date: 2019-04-09 09:29:47.244358

"""
from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision = '4a39e2cdc26a'
down_revision = 'e6ae69e9dcb9'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.create_table('experiment',
                    sa.Column('id', sa.Integer(), nullable=False),
                    sa.Column('name', sa.String(), nullable=True),
                    sa.Column('description', sa.Text(), nullable=True),
                    sa.Column('root_directory', sa.String(), nullable=True),
                    sa.Column('owner_id', sa.Integer(), nullable=True),
                    sa.Column('created_at', sa.DateTime(), nullable=True),
                    sa.ForeignKeyConstraint(['owner_id'], ['user.id'], ),
                    sa.PrimaryKeyConstraint('id')
                    )
    op.create_index(op.f('ix_experiment_id'), 'experiment', ['id'], unique=False)
    op.create_index(op.f('ix_experiment_name'), 'experiment', ['name'], unique=False)
    op.create_index(op.f('ix_experiment_owner_id'), 'experiment', ['owner_id'], unique=False)

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_index(op.f('ix_experiment_owner_id'), table_name='experiment')
    op.drop_index(op.f('ix_experiment_name'), table_name='experiment')
    op.drop_index(op.f('ix_experiment_id'), table_name='experiment')
    op.drop_table('experiment')
    # ### end Alembic commands ###
