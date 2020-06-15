<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="text-h5 primary--text">Create User</div>
      </v-card-title>
      <v-card-text>
        <template>
          <v-form v-model="valid" ref="form" lazy-validation>
            <v-text-field label="Name" v-model="name"></v-text-field>
            <v-text-field label="E-mail" type="email" v-model="email" :rules="emailRules"></v-text-field>
            <v-checkbox label="Is Admin" v-model="isAdmin"></v-checkbox>
            <v-checkbox label="Is Active" v-model="isActive"></v-checkbox>
            <v-row align="center">
              <v-col>
                <v-text-field
                  type="password"
                  ref="password"
                  label="Set Password"
                  :rules="password1Rules"
                  v-model="password1"
                ></v-text-field>
                <v-text-field
                  type="password"
                  label="Confirm Password"
                  :rules="password2Rules"
                  v-model="password2"
                ></v-text-field>
              </v-col>
            </v-row>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn @click="submit" :disabled="!valid">
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
import { userModule } from "@/modules/user";
import { IUserProfileCreate } from "@/modules/user/models";
import { required, email } from "@/utils/validators";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class CreateUser extends Vue {
  readonly userContext = userModule.context(this.$store);

  readonly emailRules = [required, email];
  readonly password1Rules = [required];
  readonly password2Rules = [required, this.passwordIsEqual];

  passwordIsEqual(v) {
    return v === this.password1 || "Password should be the same";
  }

  valid = false;
  name: string = "";
  email: string = "";
  isActive: boolean = true;
  isAdmin: boolean = false;
  password1: string = "";
  password2: string = "";

  async mounted() {
    await this.userContext.actions.getUsers();
  }

  reset() {
    this.password1 = "";
    this.password2 = "";
    this.name = "";
    this.email = "";
    this.isActive = true;
    this.isAdmin = false;
    (this.$refs.form as any).resetValidation();
  }

  cancel() {
    this.$router.back();
  }

  async submit() {
    if ((this.$refs.form as any).validate()) {
      const updatedProfile: IUserProfileCreate = {
        email: this.email,
      };
      if (this.name) {
        updatedProfile.name = this.name;
      }
      if (this.email) {
        updatedProfile.email = this.email;
      }
      updatedProfile.is_active = this.isActive;
      updatedProfile.is_admin = this.isAdmin;
      updatedProfile.password = this.password1;
      await this.userContext.actions.createUser(updatedProfile);
      this.$router.push("/main/admin/users");
    }
  }
}
</script>
