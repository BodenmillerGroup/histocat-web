"""Add gate table

Revision ID: 5a34bc7f5ed8
Revises: c18a654cf9b3
Create Date: 2020-08-14 13:23:27.886176

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY


# revision identifiers, used by Alembic.
revision = '5a34bc7f5ed8'
down_revision = 'c18a654cf9b3'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'gate',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('dataset_id', sa.Integer(), sa.ForeignKey("dataset.id", ondelete="CASCADE"), index=True, nullable=False),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('acquisition_ids', ARRAY(sa.Integer())),
        sa.Column('indices', ARRAY(sa.Integer())),
        sa.Column('cell_ids', ARRAY(sa.Integer())),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('gate')
