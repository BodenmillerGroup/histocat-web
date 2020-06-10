"""Add preset table

Revision ID: fa29b07fd05f
Revises: 35b87e175854
Create Date: 2020-06-10 17:33:04.292113

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB


# revision identifiers, used by Alembic.
revision = 'fa29b07fd05f'
down_revision = '35b87e175854'
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        'preset',
        sa.Column('id', sa.Integer(), primary_key=True, index=True),
        sa.Column('experiment_id', sa.Integer(), sa.ForeignKey("experiment.id", ondelete="CASCADE"), nullable=False, index=True),
        sa.Column('name', sa.String()),
        sa.Column('description', sa.String()),
        sa.Column('data', JSONB()),
        sa.Column('created_at', sa.DateTime(), default=sa.sql.func.now(), nullable=False),
    )


def downgrade():
    op.drop_table('preset')
